use std::io::{self, BufRead, Write};

use gridgraph_rs::{GridGraphParams, generate_dimacs};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let stdin = io::stdin();
    let mut line = String::new();
    stdin.lock().read_line(&mut line)?;

    if line.trim().is_empty() {
        return Err("expected five integers on stdin".into());
    }

    let nums: Vec<i32> = line
        .split_whitespace()
        .map(|s| s.parse::<i32>())
        .collect::<Result<Vec<_>, _>>()?;
    if nums.len() != 5 {
        return Err("expected exactly five integers".into());
    }

    let params = GridGraphParams::new(nums[0], nums[1], nums[2], nums[3], nums[4])?;
    let output = generate_dimacs(params);
    io::stdout().write_all(output.as_bytes())?;
    Ok(())
}
