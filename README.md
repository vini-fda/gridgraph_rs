# gridgraph-rs

Zero-dependency Rust rewrite of the **GRIDGRAPH** min-cost flow problem
generator. The [original](http://archive.dimacs.rutgers.edu/pub/netflow/generators/network/gridgraph/)
is a 1991 Fortran program by Mauricio G.C. Resende (AT&T Bell Labs).

The Rust binary produces **bit-for-bit identical output** to the patched
Fortran reference.

## Command-line usage

```
echo "H W MAXCAP MAXCOST SEED" | cargo run --release
```

| Parameter | Description |
|-----------|-------------|
| `H` | Grid height (rows) |
| `W` | Grid width (columns, >= 2) |
| `MAXCAP` | Maximum arc capacity |
| `MAXCOST` | Maximum arc cost |
| `SEED` | Positive integer RNG seed |

Example:

```
$ echo "3 3 100 10 12345" | cargo run --release -q
c    max-flow min-cost on a  3 by   3 grid
c    cap: [0,     100]   cost: [0,      10]   seed:     12345
p min              11             18
n              10             40
n              11            -40
a          1         2 0         10         8
...
```

## Graph structure

- `H * W + 2` nodes: an `H x W` grid, a source, and a sink.
- `2 * (H-1) * W + W + H` directed arcs: right and down within the grid,
  source to first column, last column to sink.
- Capacities and costs drawn uniformly from `[1, MAXCAP]` and `[1, MAXCOST]`.
- Output is DIMACS `.min` format with supply/demand set to the max s-t flow
  (computed via Dinic's algorithm).

## Repository structure

```
src/main.rs              Rust implementation
gridgraph_original/
  ggraph1.f              Original Fortran source (unmodified)
  Makefile               Builds the Fortran reference binary
```

## Building & testing

```
cargo build --release
cargo test               # unit + end-to-end (compares against Fortran reference)
```

The end-to-end tests build the Fortran binary automatically via `make`
(requires `gfortran` and `perl`).

## References

- M.G.C. Resende, "GRIDGRAPH generator", AT&T Bell Labs, 1991.
- L. Schrage, "A More Portable FORTRAN Random Number Generator", *ACM TOMS*, 1979.
- D. Goldfarb & M.D. Grigoriadis, "A computational comparison of the Dinic
  and network simplex methods for maximum flow", *Annals of Operations Research* 13, 1988.

## Library API

`gridgraph_rs` is now a reusable library crate. Add it to another project's
`Cargo.toml` with

```
cargo add gridgraph_rs
```

```rust
use gridgraph_rs::{generate_instance, GridGraphParams};

fn main() {
    let params = GridGraphParams::new(3, 3, 100, 10, 12345).unwrap();
    let instance = generate_instance(params);

    // Render DIMACS output
    println!("{}", instance.to_dimacs_string());

    // Or inspect individual arcs
    for arc in instance.arcs() {
        println!("{} -> {} (cap {}, cost {})", arc.from, arc.to, arc.capacity, arc.cost);
    }
}
```

Public API highlights:

- `GridGraphParams::new` validates the five GRIDGRAPH parameters.
- `generate_instance` returns a [`DimacsInstance`](src/lib.rs) with metadata and arc list.
- `DimacsInstance::write_dimacs` / `to_dimacs_string` format DIMACS output.
- `generate_dimacs` is a convenience helper returning the formatted string in one step.
