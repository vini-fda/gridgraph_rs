#!/usr/bin/env -S just --justfile

set shell := ["bash", "-euo", "pipefail", "-c"]

# Core workflows
fmt:
    cargo fmt --all

fmt-check:
    cargo fmt --all -- --check

clippy:
    cargo clippy --workspace --all-targets -- -D warnings

test:
    cargo test

build profile='release':
    cargo build --workspace --profile {{ profile }}

clean:
    cargo clean

ci: fmt-check clippy test
