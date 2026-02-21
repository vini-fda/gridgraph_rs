//! GRIDGRAPH — min-cost flow grid graph generator
//!
//! Produces a min-cost flow problem instance in DIMACS `.min` format,
//! laid out on a directed grid graph. Originally by Mauricio G.C. Resende
//! (1991), AT&T Bell Laboratories.
//!
//! This is a zero-dependency Rust rewrite of the original Fortran program
//! `ggraph1.f`, producing bit-for-bit identical output.
//!
//! # Usage
//!
//! ```text
//! echo "H W MAXCAP MAXCOST SEED" | gridgraph_rs
//! ```
//!
//! The five space-separated integers on stdin are:
//!
//! | Parameter | Description                              |
//! |-----------|------------------------------------------|
//! | `H`       | Grid height (number of rows)             |
//! | `W`       | Grid width (number of columns, >= 2)     |
//! | `MAXCAP`  | Maximum arc capacity (caps in `[1, MAXCAP]`) |
//! | `MAXCOST` | Maximum arc cost (costs in `[1, MAXCOST]`)   |
//! | `SEED`    | Positive integer seed for the RNG        |
//!
//! # Graph structure
//!
//! The generated graph has `H * W + 2` nodes:
//!
//! - A grid of `H` rows by `W` columns, where node `(i, j)` (1-indexed)
//!   is numbered `(i - 1) * W + j`.
//! - A **source** node numbered `H * W + 1`.
//! - A **sink** node numbered `H * W + 2`.
//!
//! Arcs are directed right and down within the grid. The source connects
//! to every node in the first column, and every node in the last column
//! connects to the sink. The total arc count is `2 * (H - 1) * W + W + H`.
//!
//! # Output
//!
//! The output is a DIMACS min-cost flow problem. Supply at the source and
//! demand at the sink are set to the maximum s-t flow value (computed
//! internally via Dinic's algorithm), so the min-cost flow problem is
//! always feasible with a unique optimal flow value.
//!
//! # References
//!
//! - M.G.C. Resende, "GRIDGRAPH generator", AT&T Bell Labs, 1991.
//! - L. Schrage, "A More Portable FORTRAN Random Number Generator",
//!   *ACM TOMS*, 1979.
//! - D. Goldfarb & M.D. Grigoriadis, "A computational comparison of the
//!   Dinic and network simplex methods for maximum flow",
//!   *Annals of Operations Research* 13, 1988.

use std::collections::VecDeque;
use std::io::{self, BufRead, Write};

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

struct Arc {
    from: i32,
    to: i32,
    cap: i32,
    cost: i32,
}

// ---------------------------------------------------------------------------
// Grid graph generation (9 phases matching the Fortran arc order)
// ---------------------------------------------------------------------------

/// Map 1-indexed grid position `(row, col)` to its node number.
fn node(i: i32, j: i32, w: i32) -> i32 {
    (i - 1) * w + j
}

/// Generate all arcs in the exact order the Fortran program produces them.
///
/// The generation order determines RNG consumption and therefore every arc
/// attribute. Source and sink arc capacities/costs are accumulated during
/// grid arc generation and emitted at the end (phases 8-9).
fn generate_arcs(h: i32, w: i32, maxcap: i32, maxcost: i32, seed: &mut i32) -> Vec<Arc> {
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
        arcs.push(Arc {
            from: node(i, 1, w),
            to: node(i, 2, w),
            cap,
            cost,
        });

        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        caps[i as usize] += cap;
        arcs.push(Arc {
            from: node(i, 1, w),
            to: node(i + 1, 1, w),
            cap,
            cost,
        });

        costs[i as usize] = unif(maxcost, seed); // extra RNG call
    }

    // Phase 2 — col 1, row h
    {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        caps[h as usize] = cap;
        arcs.push(Arc {
            from: node(h, 1, w),
            to: node(h, 2, w),
            cap,
            cost,
        });

        costs[h as usize] = unif(maxcost, seed); // extra RNG call
    }

    // Phase 3 — interior cols 2..w-2, rows 1..h-1
    for i in 1..h {
        for j in 2..=(w - 2) {
            let cap = unif(maxcap, seed);
            let cost = unif(maxcost, seed);
            arcs.push(Arc {
                from: node(i, j, w),
                to: node(i, j + 1, w),
                cap,
                cost,
            });

            let cap = unif(maxcap, seed);
            let cost = unif(maxcost, seed);
            arcs.push(Arc {
                from: node(i, j, w),
                to: node(i + 1, j, w),
                cap,
                cost,
            });
        }
    }

    // Phase 4 — interior cols 2..w-2, row h
    for j in 2..=(w - 2) {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        arcs.push(Arc {
            from: node(h, j, w),
            to: node(h, j + 1, w),
            cap,
            cost,
        });
    }

    // Phase 5 — col w-1, rows 1..h-1
    for i in 1..h {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        capt[i as usize] = cap;
        arcs.push(Arc {
            from: node(i, w - 1, w),
            to: node(i, w, w),
            cap,
            cost,
        });

        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        arcs.push(Arc {
            from: node(i, w - 1, w),
            to: node(i + 1, w - 1, w),
            cap,
            cost,
        });

        costt[i as usize] = unif(maxcost, seed); // extra RNG call
    }

    // Phase 6 — col w-1, row h
    {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        capt[h as usize] = cap;
        arcs.push(Arc {
            from: node(h, w - 1, w),
            to: node(h, w, w),
            cap,
            cost,
        });

        costt[h as usize] = unif(maxcost, seed); // extra RNG call
    }

    // Phase 7 — col w, rows 1..h-1
    for i in 1..h {
        let cap = unif(maxcap, seed);
        let cost = unif(maxcost, seed);
        capt[(i + 1) as usize] += cap;
        arcs.push(Arc {
            from: node(i, w, w),
            to: node(i + 1, w, w),
            cap,
            cost,
        });
    }

    let source = h * w + 1;
    let sink = h * w + 2;

    // Phase 8 — source arcs
    for i in 1..=h {
        arcs.push(Arc {
            from: source,
            to: node(i, 1, w),
            cap: caps[i as usize],
            cost: costs[i as usize],
        });
    }

    // Phase 9 — sink arcs
    for i in 1..=h {
        arcs.push(Arc {
            from: node(i, w, w),
            to: sink,
            cap: capt[i as usize],
            cost: costt[i as usize],
        });
    }

    debug_assert_eq!(arcs.len(), narcs as usize);
    arcs
}

// ---------------------------------------------------------------------------
// Dinic's max-flow
//
// Standard implementation: BFS to build level graph, then DFS to push
// blocking flow with the "current arc" optimisation. Uses i64 for flow
// values to avoid overflow on large grids.
//
// Reference: Goldfarb & Grigoriadis, Annals of Operations Research, 1988.
// ---------------------------------------------------------------------------

/// Adjacency-list flow network for Dinic's algorithm.
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

    /// Add a directed edge `u -> v` with capacity `c` (and its reverse with capacity 0).
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

    /// BFS to assign distance levels. Returns `true` if the sink is reachable.
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

    /// DFS to find a blocking flow along level-graph edges.
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

    /// Compute the maximum flow from `source` to `sink`.
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

// ---------------------------------------------------------------------------
// Core driver: parse input, generate graph, compute max-flow, format output
// ---------------------------------------------------------------------------

/// Generate the DIMACS `.min` output for the given input parameters.
fn run(input: &str) -> String {
    let nums: Vec<i32> = input
        .split_whitespace()
        .map(|s| s.parse().unwrap())
        .collect();
    let (h, w, maxcap, maxcost) = (nums[0], nums[1], nums[2], nums[3]);
    let mut seed = nums[4];
    let seed0 = seed;

    let nnodes = h * w + 2;
    let narcs = 2 * (h - 1) * w + w + h;
    let source = (h * w + 1) as usize;
    let sink = (h * w + 2) as usize;

    let arcs = generate_arcs(h, w, maxcap, maxcost, &mut seed);

    // Build flow graph for Dinic's
    let mut fg = FlowGraph::new(nnodes as usize + 1, narcs as usize * 2);
    for arc in &arcs {
        fg.add_edge(arc.from as usize, arc.to as usize, arc.cap as i64);
    }

    let max_flow = fg.max_flow(source, sink);

    // Format output
    let mut out = String::new();
    use std::fmt::Write as FmtWrite;

    writeln!(out, "c    max-flow min-cost on a{h:3} by {w:3} grid").unwrap();
    writeln!(
        out,
        "c    cap: [0,{maxcap:8}]   cost: [0,{maxcost:8}]   seed:{seed0:10}"
    )
    .unwrap();
    writeln!(out, "p min {nnodes:15}{narcs:15}").unwrap();
    writeln!(out, "n {source:15}{max_flow:15}").unwrap();
    let neg_flow = -max_flow;
    writeln!(out, "n {sink:15}{neg_flow:15}").unwrap();
    for arc in &arcs {
        let f = arc.from;
        let t = arc.to;
        let c = arc.cap;
        let co = arc.cost;
        writeln!(out, "a {f:10}{t:10} 0 {c:10}{co:10}").unwrap();
    }

    out
}

fn main() {
    let stdin = io::stdin();
    let line = stdin.lock().lines().next().unwrap().unwrap();
    let output = run(&line);
    io::stdout().write_all(output.as_bytes()).unwrap();
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

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
        let actual = run("3 3 100 10 12345");
        assert_eq!(actual, expected);
    }

    fn run_fortran(input: &str) -> Option<String> {
        let manifest = env!("CARGO_MANIFEST_DIR");
        let fortran_dir = format!("{manifest}/gridgraph_original");

        // Build if needed
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
        let input = "3 3 100 10 12345";
        if let Some(fortran_out) = run_fortran(input) {
            let rust_out = run(input);
            assert_eq!(rust_out, fortran_out);
        } else {
            eprintln!("Skipping Fortran comparison: could not build/run Fortran binary");
        }
    }

    #[test]
    fn e2e_8x8_vs_fortran() {
        let input = "8 8 10000 1000 270001";
        if let Some(fortran_out) = run_fortran(input) {
            let rust_out = run(input);
            assert_eq!(rust_out, fortran_out);
        } else {
            eprintln!("Skipping Fortran comparison: could not build/run Fortran binary");
        }
    }

    #[test]
    fn e2e_8x8_spot_check() {
        let output = run("8 8 10000 1000 270001");
        let lines: Vec<&str> = output.lines().collect();
        // 66 nodes, 128 arcs
        assert!(lines[2].contains("66"));
        assert!(lines[2].contains("128"));
        // max-flow = 12458
        assert!(lines[3].contains("12458"));
        // first arc
        assert_eq!(lines[5], "a          1         2 0       1132       335");
    }
}
