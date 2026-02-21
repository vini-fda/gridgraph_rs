use std::io::{self, BufRead, BufWriter, Write};

use gridgraph_rs::{GridGraphParams, generate_instance};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let nums = if args.is_empty() {
        read_numbers_from_stdin()?
    } else {
        let nums = parse_numbers(args.iter().map(|s| s.as_str()))?;
        if nums.len() != 5 {
            return Err("expected exactly five integers".into());
        }
        nums
    };

    let params = GridGraphParams::new(nums[0], nums[1], nums[2], nums[3], nums[4])?;
    let instance = generate_instance(params);
    let stdout = io::stdout();
    let mut writer = BufWriter::new(stdout.lock());
    instance.write_dimacs_io(&mut writer)?;
    writer.flush()?;
    Ok(())
}

fn read_numbers_from_stdin() -> Result<Vec<i32>, Box<dyn std::error::Error>> {
    let stdin = io::stdin();
    let mut line = String::new();
    stdin.lock().read_line(&mut line)?;
    if line.trim().is_empty() {
        return Err("expected five integers on stdin".into());
    }
    let nums = parse_numbers(line.split_whitespace())?;
    if nums.len() != 5 {
        return Err("expected exactly five integers".into());
    }
    Ok(nums)
}

fn parse_numbers<I, S>(inputs: I) -> Result<Vec<i32>, std::num::ParseIntError>
where
    I: IntoIterator<Item = S>,
    S: AsRef<str>,
{
    inputs
        .into_iter()
        .map(|s| s.as_ref().parse::<i32>())
        .collect()
}
