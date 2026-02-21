//! GRIDGRAPH, a min-cost flow grid graph generator
//!
//! Produces a min-cost flow problem instance in DIMACS `.min` format,
//! laid out on a directed grid graph. Originally by Mauricio G.C. Resende
//! (1991), AT&T Bell Laboratories.
//!
//! This is a Rust rewrite of the [original](http://archive.dimacs.rutgers.edu/pub/netflow/generators/network/gridgraph/)
//! Fortran program `ggraph1.f`.
//!
//! # Usage
//!
//! ```text
//! use gridgraph_rs::{generate_instance, GridGraphParams};
//!
//! let params = GridGraphParams::new(3, 3, 100, 10, 12345).unwrap();
//! let instance = generate_instance(params);
//! println!("{}", instance.to_dimacs_string());
//! ```
//!
//! The five integers correspond to grid dimensions, capacity/cost bounds
//! and RNG seed. See [`GridGraphParams`] for details.

use std::collections::VecDeque;
use std::fmt;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Parameters that fully describe a GRIDGRAPH instance.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct GridGraphParams {
    /// Grid height (rows), must be positive.
    pub height: i32,
    /// Grid width (columns), must be at least 2.
    pub width: i32,
    /// Maximum arc capacity, sampled uniformly from `[1, max_capacity]`.
    pub max_capacity: i32,
    /// Maximum arc cost, sampled uniformly from `[1, max_cost]`.
    pub max_cost: i32,
    /// Positive integer seed for the Park-Miller RNG.
    pub seed: i32,
}

impl GridGraphParams {
    /// Validate and build a new parameter set.
    pub fn new(
        height: i32,
        width: i32,
        max_capacity: i32,
        max_cost: i32,
        seed: i32,
    ) -> Result<Self, GridGraphError> {
        if height <= 0 {
            return Err(GridGraphError::InvalidHeight);
        }
        if width < 2 {
            return Err(GridGraphError::InvalidWidth);
        }
        if max_capacity <= 0 {
            return Err(GridGraphError::InvalidMaxCapacity);
        }
        if max_cost <= 0 {
            return Err(GridGraphError::InvalidMaxCost);
        }
        if seed <= 0 {
            return Err(GridGraphError::InvalidSeed);
        }
        Ok(Self {
            height,
            width,
            max_capacity,
            max_cost,
            seed,
        })
    }
}

/// Errors returned when the supplied parameters are invalid.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GridGraphError {
    InvalidHeight,
    InvalidWidth,
    InvalidMaxCapacity,
    InvalidMaxCost,
    InvalidSeed,
}

impl fmt::Display for GridGraphError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use GridGraphError::*;
        match self {
            InvalidHeight => write!(f, "height must be positive"),
            InvalidWidth => write!(f, "width must be >= 2"),
            InvalidMaxCapacity => write!(f, "max_capacity must be positive"),
            InvalidMaxCost => write!(f, "max_cost must be positive"),
            InvalidSeed => write!(f, "seed must be positive"),
        }
    }
}

impl std::error::Error for GridGraphError {}

/// Directed arc within the generated grid graph.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GridArc {
    /// Tail node identifier (1-indexed for grid nodes).
    pub from: i32,
    /// Head node identifier.
    pub to: i32,
    /// Arc capacity.
    pub capacity: i32,
    /// Arc cost.
    pub cost: i32,
}

/// Complete DIMACS min-cost flow instance for a GRIDGRAPH.
#[derive(Debug, Clone)]
pub struct DimacsInstance {
    params: GridGraphParams,
    arcs: Vec<GridArc>,
    node_count: i32,
    arc_count: i32,
    source: i32,
    sink: i32,
    max_flow: i64,
}

impl DimacsInstance {
    /// Original parameters used to generate this instance.
    pub fn params(&self) -> GridGraphParams {
        self.params
    }

    /// Total number of nodes (`H * W + 2`).
    pub fn node_count(&self) -> i32 {
        self.node_count
    }

    /// Total number of arcs.
    pub fn arc_count(&self) -> i32 {
        self.arc_count
    }

    /// Source node identifier (`H * W + 1`).
    pub fn source(&self) -> i32 {
        self.source
    }

    /// Sink node identifier (`H * W + 2`).
    pub fn sink(&self) -> i32 {
        self.sink
    }

    /// Maximum s-t flow value used for supply/demand.
    pub fn max_flow(&self) -> i64 {
        self.max_flow
    }

    /// Borrow the arc list in Fortran order.
    pub fn arcs(&self) -> &[GridArc] {
        &self.arcs
    }

    /// Render the DIMACS `.min` output to a writer.
    pub fn write_dimacs<W: fmt::Write>(&self, mut writer: W) -> fmt::Result {
        let params = self.params;
        writeln!(
            writer,
            "c    max-flow min-cost on a{h:3} by {w:3} grid",
            h = params.height,
            w = params.width
        )?;
        writeln!(
            writer,
            "c    cap: [0,{cap:8}]   cost: [0,{cost:8}]   seed:{seed:10}",
            cap = params.max_capacity,
            cost = params.max_cost,
            seed = params.seed
        )?;
        writeln!(
            writer,
            "p min {nodes:15}{arcs:15}",
            nodes = self.node_count,
            arcs = self.arc_count
        )?;
        writeln!(
            writer,
            "n {src:15}{flow:15}",
            src = self.source,
            flow = self.max_flow
        )?;
        writeln!(
            writer,
            "n {sink:15}{neg_flow:15}",
            sink = self.sink,
            neg_flow = -self.max_flow
        )?;
        for arc in &self.arcs {
            writeln!(
                writer,
                "a {from:10}{to:10} 0 {cap:10}{cost:10}",
                from = arc.from,
                to = arc.to,
                cap = arc.capacity,
                cost = arc.cost
            )?;
        }
        Ok(())
    }

    /// Convenience method returning the formatted DIMACS string.
    pub fn to_dimacs_string(&self) -> String {
        let mut out = String::new();
        self.write_dimacs(&mut out).unwrap();
        out
    }
}

/// Generate a grid graph instance along with all derived metadata.
pub fn generate_instance(params: GridGraphParams) -> DimacsInstance {
    let h = params.height;
    let w = params.width;
    let maxcap = params.max_capacity;
    let maxcost = params.max_cost;

    let nnodes = h * w + 2;
    let narcs = 2 * (h - 1) * w + w + h;
    let source = h * w + 1;
    let sink = h * w + 2;

    let mut seed = params.seed;
    let arcs = generate_arcs(h, w, maxcap, maxcost, &mut seed);

    // Build flow graph for Dinic's algorithm
    let mut fg = FlowGraph::new(nnodes as usize + 1, narcs as usize * 2);
    for arc in &arcs {
        fg.add_edge(arc.from as usize, arc.to as usize, arc.capacity as i64);
    }
    let max_flow = fg.max_flow(source as usize, sink as usize);

    DimacsInstance {
        params,
        arcs,
        node_count: nnodes,
        arc_count: narcs,
        source,
        sink,
        max_flow,
    }
}

/// Generate a DIMACS `.min` string for the supplied parameters.
pub fn generate_dimacs(params: GridGraphParams) -> String {
    generate_instance(params).to_dimacs_string()
}

// ---------------------------------------------------------------------------
// RNG — Schrage's Park-Miller LCG
//
// Uses f32 arithmetic throughout to match Fortran's 32-bit REAL.
// Reference: L. Schrage, ACM TOMS, 1979.
// ---------------------------------------------------------------------------

const A: i32 = 16807;
const P: i32 = 2_147_483_647; // 2^31 - 1
const B15: i32 = 32768; // 2^15
const B16: i32 = 65536; // 2^16

/// Advance the Park-Miller LCG state and return a value in [0, 1).
fn schrage_rand(seed: &mut i32) -> f32 {
    let ix = *seed;
    let xhi = ix / B16;
    let xalo = (ix - xhi * B16) * A;
    let leftlo = xalo / B16;
    let fhi = xhi * A + leftlo;
    let k = fhi / B15;
    let mut new_ix = ((xalo - leftlo * B16) - P) + (fhi - k * B15) * B16 + k;
    if new_ix < 0 {
        new_ix += P;
    }
    *seed = new_ix;
    #[allow(clippy::excessive_precision)]
    let scale = 4.656_612_875e-10_f32;
    new_ix as f32 * scale
}

/// Return a uniform random integer in `[1, max]`.
fn unif(max: i32, seed: &mut i32) -> i32 {
    (1.0_f32 + schrage_rand(seed) * (max - 1) as f32) as i32
}

// ---------------------------------------------------------------------------
// Arc representation
// ---------------------------------------------------------------------------

fn node(i: i32, j: i32, w: i32) -> i32 {
    (i - 1) * w + j
}

fn generate_arcs(h: i32, w: i32, maxcap: i32, maxcost: i32, seed: &mut i32) -> Vec<GridArc> {
    let narcs = 2 * (h - 1) * w + w + h;
    let mut arcs = Vec::with_capacity(narcs as usize);

    let mut caps = vec![0i32; (h + 1) as usize]; // 1-indexed
    let mut costs = vec![0i32; (h + 1) as usize];
    let mut capt = vec![0i32; (h + 1) as usize];
    let mut costt = vec![0i32; (h + 1) as usize];

    // Phase 1 — col 1, rows 1..h-1
    for i in 1..h {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        caps[i as usize] = cap;
        arcs.push(GridArc {
            from: node(i, 1, w),
            to: node(i, 2, w),
            capacity: cap,
            cost,
        });

        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        caps[i as usize] += cap;
        arcs.push(GridArc {
            from: node(i, 1, w),
            to: node(i + 1, 1, w),
            capacity: cap,
            cost,
        });

        costs[i as usize] = unif(maxcost, seed); // extra RNG call
    }

    // Phase 2 — col 1, row h
    {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        caps[h as usize] = cap;
        arcs.push(GridArc {
            from: node(h, 1, w),
            to: node(h, 2, w),
            capacity: cap,
            cost,
        });

        costs[h as usize] = unif(maxcost, seed); // extra RNG call
    }

    // Phase 3 — interior cols 2..w-2, rows 1..h-1
    for i in 1..h {
        for j in 2..=(w - 2) {
            let cap = unif(maxcap, seed);
            let cost = unif(maxcost, seed);
            arcs.push(GridArc {
                from: node(i, j, w),
                to: node(i, j + 1, w),
                capacity: cap,
                cost,
            });

            let cap = unif(maxcap, seed);
            let cost = unif(maxcost, seed);
            arcs.push(GridArc {
                from: node(i, j, w),
                to: node(i + 1, j, w),
                capacity: cap,
                cost,
            });
        }
    }

    // Phase 4 — interior cols 2..w-2, row h
    for j in 2..=(w - 2) {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        arcs.push(GridArc {
            from: node(h, j, w),
            to: node(h, j + 1, w),
            capacity: cap,
            cost,
        });
    }

    // Phase 5 — col w-1, rows 1..h-1
    for i in 1..h {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        capt[i as usize] = cap;
        arcs.push(GridArc {
            from: node(i, w - 1, w),
            to: node(i, w, w),
            capacity: cap,
            cost,
        });

        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        arcs.push(GridArc {
            from: node(i, w - 1, w),
            to: node(i + 1, w - 1, w),
            capacity: cap,
            cost,
        });

        costt[i as usize] = unif(maxcost, seed); // extra RNG call
    }

    // Phase 6 — col w-1, row h
    {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        capt[h as usize] = cap;
        arcs.push(GridArc {
            from: node(h, w - 1, w),
            to: node(h, w, w),
            capacity: cap,
            cost,
        });

        costt[h as usize] = unif(maxcost, seed); // extra RNG call
    }

    // Phase 7 — col w, rows 1..h-1
    for i in 1..h {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        capt[(i + 1) as usize] += cap;
        arcs.push(GridArc {
            from: node(i, w, w),
            to: node(i + 1, w, w),
            capacity: cap,
            cost,
        });
    }

    let source = h * w + 1;
    let sink = h * w + 2;

    // Phase 8 — source arcs
    for i in 1..=h {
        arcs.push(GridArc {
            from: source,
            to: node(i, 1, w),
            capacity: caps[i as usize],
            cost: costs[i as usize],
        });
    }

    // Phase 9 — sink arcs
    for i in 1..=h {
        arcs.push(GridArc {
            from: node(i, w, w),
            to: sink,
            capacity: capt[i as usize],
            cost: costt[i as usize],
        });
    }

    debug_assert_eq!(arcs.len(), narcs as usize);
    arcs
}

// ---------------------------------------------------------------------------
// Dinic's max-flow
// ---------------------------------------------------------------------------

struct FlowGraph {
    n: usize,
    head: Vec<i32>,
    to: Vec<i32>,
    cap: Vec<i64>,
    next: Vec<i32>,
}

impl FlowGraph {
    fn new(n: usize, arc_hint: usize) -> Self {
        Self {
            n,
            head: vec![-1; n],
            to: Vec::with_capacity(arc_hint),
            cap: Vec::with_capacity(arc_hint),
            next: Vec::with_capacity(arc_hint),
        }
    }

    fn add_edge(&mut self, u: usize, v: usize, c: i64) {
        let idx = self.to.len() as i32;
        self.to.push(v as i32);
        self.cap.push(c);
        self.next.push(self.head[u]);
        self.head[u] = idx;

        let idx = self.to.len() as i32;
        self.to.push(u as i32);
        self.cap.push(0);
        self.next.push(self.head[v]);
        self.head[v] = idx;
    }

    fn bfs(&self, source: usize, sink: usize, level: &mut [i32]) -> bool {
        level.iter_mut().for_each(|l| *l = -1);
        level[source] = 0;
        let mut queue = VecDeque::new();
        queue.push_back(source);
        while let Some(u) = queue.pop_front() {
            let mut e = self.head[u];
            while e != -1 {
                let v = self.to[e as usize] as usize;
                if level[v] == -1 && self.cap[e as usize] > 0 {
                    level[v] = level[u] + 1;
                    queue.push_back(v);
                }
                e = self.next[e as usize];
            }
        }
        level[sink] != -1
    }

    fn dfs(&mut self, u: usize, sink: usize, pushed: i64, level: &[i32], iter: &mut [i32]) -> i64 {
        if u == sink {
            return pushed;
        }
        while iter[u] != -1 {
            let e = iter[u] as usize;
            let v = self.to[e] as usize;
            if level[v] == level[u] + 1 && self.cap[e] > 0 {
                let d = self.dfs(v, sink, pushed.min(self.cap[e]), level, iter);
                if d > 0 {
                    self.cap[e] -= d;
                    self.cap[e ^ 1] += d;
                    return d;
                }
            }
            iter[u] = self.next[iter[u] as usize];
        }
        0
    }

    fn max_flow(&mut self, source: usize, sink: usize) -> i64 {
        let mut flow = 0i64;
        let mut level = vec![0i32; self.n];
        let mut iter = vec![0i32; self.n];
        while self.bfs(source, sink, &mut level) {
            iter.copy_from_slice(&self.head);
            loop {
                let f = self.dfs(source, sink, i64::MAX, &level, &mut iter);
                if f == 0 {
                    break;
                }
                flow += f;
            }
        }
        flow
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::process::Command;

    #[test]
    fn test_schrage_rand_seed_270001() {
        let mut seed: i32 = 270_001;
        let r1 = schrage_rand(&mut seed);
        assert_eq!(seed, 242_939_513);
        assert!((r1 - 0.11313).abs() < 1e-4, "r1 = {r1}");

        let r2 = schrage_rand(&mut seed);
        assert_eq!(seed, 717_982_044);
        assert!((r2 - 0.33434).abs() < 1e-4, "r2 = {r2}");
    }

    #[test]
    fn test_unif() {
        let mut seed: i32 = 270_001;
        assert_eq!(unif(10_000, &mut seed), 1132);
        assert_eq!(unif(1_000, &mut seed), 335);
    }

    #[test]
    fn test_arc_count() {
        assert_eq!(2 * (3 - 1) * 3 + 3 + 3, 18);
        assert_eq!(2 * (8 - 1) * 8 + 8 + 8, 128);
    }

    #[test]
    fn test_node_numbering() {
        assert_eq!(node(1, 1, 3), 1);
        assert_eq!(node(2, 3, 3), 6);
        assert_eq!(node(3, 3, 3), 9);
    }

    #[test]
    fn e2e_3x3_known_output() {
        let params = GridGraphParams::new(3, 3, 100, 10, 12345).unwrap();
        let expected = "\
c    max-flow min-cost on a  3 by   3 grid
c    cap: [0,     100]   cost: [0,      10]   seed:     12345
p min              11             18
n              10             40
n              11            -40
a          1         2 0         10         8
a          1         4 0         94         1
a          4         5 0          6         7
a          4         7 0         58         9
a          7         8 0         33         2
a          2         3 0         79         9
a          2         5 0         15         4
a          5         6 0         27         2
a          5         8 0         86         2
a          8         9 0         24         9
a          3         6 0         31         6
a          6         9 0         51         2
a         10         1 0        104         1
a         10         4 0         64         8
a         10         7 0         33         3
a          3        11 0         79         1
a          6        11 0         58         9
a          9        11 0         75         9
";
        let actual = generate_dimacs(params);
        assert_eq!(actual, expected);
    }

    fn run_fortran(input: &str) -> Option<String> {
        let manifest = env!("CARGO_MANIFEST_DIR");
        let fortran_dir = format!("{manifest}/gridgraph_original");

        let status = Command::new("make")
            .args(["-C", &fortran_dir])
            .status()
            .ok()?;
        if !status.success() {
            return None;
        }

        let binary = format!("{fortran_dir}/ggraph");
        let output = Command::new(&binary)
            .stdin(std::process::Stdio::piped())
            .stdout(std::process::Stdio::piped())
            .stderr(std::process::Stdio::piped())
            .spawn()
            .and_then(|mut child| {
                use std::io::Write;
                child
                    .stdin
                    .as_mut()
                    .unwrap()
                    .write_all(input.as_bytes())
                    .unwrap();
                child.wait_with_output()
            })
            .ok()?;

        if !output.status.success() {
            return None;
        }
        Some(String::from_utf8(output.stdout).unwrap())
    }

    #[test]
    fn e2e_3x3_vs_fortran() {
        let params = GridGraphParams::new(3, 3, 100, 10, 12345).unwrap();
        let input = "3 3 100 10 12345";
        if let Some(fortran_out) = run_fortran(input) {
            let rust_out = generate_dimacs(params);
            assert_eq!(rust_out, fortran_out);
        } else {
            eprintln!("Skipping Fortran comparison: could not build/run Fortran binary");
        }
    }

    #[test]
    fn e2e_8x8_vs_fortran() {
        let params = GridGraphParams::new(8, 8, 10_000, 1_000, 270_001).unwrap();
        let input = "8 8 10000 1000 270001";
        if let Some(fortran_out) = run_fortran(input) {
            let rust_out = generate_dimacs(params);
            assert_eq!(rust_out, fortran_out);
        } else {
            eprintln!("Skipping Fortran comparison: could not build/run Fortran binary");
        }
    }

    #[test]
    fn e2e_8x8_spot_check() {
        let params = GridGraphParams::new(8, 8, 10_000, 1_000, 270_001).unwrap();
        let output = generate_dimacs(params);
        let lines: Vec<&str> = output.lines().collect();
        assert!(lines[2].contains("66"));
        assert!(lines[2].contains("128"));
        assert!(lines[3].contains("12458"));
        assert_eq!(lines[5], "a          1         2 0       1132       335");
    }
}
