extern crate tempfile;
#[macro_use]
extern crate difference;

use std::env;
use std::fs::File;
use std::io::Read;
use std::process::Command;

fn test_resume_with(total_iters: u64, first_iters: u64) {
    let dir = tempfile::tempdir().expect("Unable to create temp directory");
    let mut root = env::current_exe()
        .unwrap()
        .parent()
        .expect("executable's directory")
        .to_path_buf();
    if root.ends_with("deps") {
        root.pop();
    }
    println!("root and dir are {:?} and {:?}", root, dir);
    let mut cmd = Command::new(root.join("histogram"));
    cmd.env("RUST_BACKTRACE", "1");
    cmd.current_dir(dir.path()).args(&["--help"]);
    let out = cmd.output().expect("command failed to run");
    println!("{}", String::from_utf8_lossy(&out.stdout));
    println!("{}", String::from_utf8_lossy(&out.stderr));
    assert!(out.status.success());

    let common_flags = &[
        "--sw-N=100",
        "--sw-filling-fraction=0.3",
        "--sw-well-width=1.3",
        "--sad-min-T=0.5",
        "--acceptance-rate=0.5",
    ];

    println!("About to start the big guy for {} moves", total_iters);
    let mut cmd = Command::new(root.join("histogram"));
    let total_max_iter = format!("--max-iter={}", total_iters);
    let first_max_iter = format!("--max-iter={}", first_iters);
    cmd.env("RUST_BACKTRACE", "1");
    cmd.current_dir(dir.path())
        .args(common_flags)
        .args(&[&total_max_iter, "--save-as=big-guy.yaml"]);
    let out = cmd.output().expect("command failed to run");
    println!("{}", String::from_utf8_lossy(&out.stdout));
    println!("{}", String::from_utf8_lossy(&out.stderr));
    assert!(out.status.success());
    println!("FINISHED BIG SIMULATION\n\n\n");

    println!("About to start the small guy ifrst");
    let mut cmd = Command::new(root.join("histogram"));
    cmd.current_dir(dir.path())
        .args(common_flags)
        .args(&[&first_max_iter, "--save-as=small-guy.yaml"]);
    let out = cmd.output().expect("command failed to run");
    println!("{}", String::from_utf8_lossy(&out.stdout));
    println!("{}", String::from_utf8_lossy(&out.stderr));
    assert!(out.status.success());
    println!("FINISHED SHORT SIMULATION\n\n\n");

    println!("About to start the small guy");
    let mut cmd = Command::new(root.join("histogram"));
    cmd.current_dir(dir.path())
        .args(common_flags)
        .args(&[&total_max_iter, "--save-as=small-guy.yaml"]);
    let out = cmd.output().expect("command failed to run");
    println!("{}", String::from_utf8_lossy(&out.stdout));
    println!("{}", String::from_utf8_lossy(&out.stderr));
    assert!(out.status.success());
    println!("FINISHED RESUMED SIMULATION\n\n\n");

    let mut f1 = File::open(dir.path().join("big-guy.yaml")).unwrap();
    let mut s1 = String::new();
    f1.read_to_string(&mut s1).unwrap();

    let mut f2 = File::open(dir.path().join("small-guy.yaml")).unwrap();
    let mut s2 = String::new();
    f2.read_to_string(&mut s2).unwrap();

    println!("\n\n\nRUNNING WITH {} and {}", total_iters, first_iters);
    assert_diff!(&s1, &s2, "\n", 2); // This corresponds to save_as differing.
}

#[test]
fn test_resume_short() {
    test_resume_with(2, 1);
    test_resume_with(5, 4);
    test_resume_with(1_000, 1);
    test_resume_with(1_000, 999);
    test_resume_with(1_000, 500);
}

#[test]
fn test_resume_long() {
    test_resume_with(1_000_000, 500_000);
}
