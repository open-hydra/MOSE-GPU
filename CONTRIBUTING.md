# Contributing to MOSE-GPU

Thanks for your interest in MOSE-GPU. This page is the short version: branch, work, open a pull request. The
[Developer Guide][devguide] has the detail.

## Workflow

Every change reaches `main` through a pull request.

```bash
# 1. Clone with submodules
git clone --recurse-submodules https://github.com/open-hydra/MOSE-GPU.git
cd MOSE-GPU

# 2. Branch from main
git switch main && git pull
git switch -c feature/my-change

# 3. Build
./install.sh build --production
./install.sh compile          # iterative rebuilds after the first build

# 4. Work, in small focused commits
git commit -m "Port the SA-RC rotation correction"

# 5. Push and open a pull request against main
git push -u origin feature/my-change
```

## Before you open the pull request

MOSE-GPU must keep giving the results of MOSE. Run the regression table on the test cases of MOSE ([Testing][testing]):

```bash
python3 test/run_mose_tests.py stage --mose /path/to/MOSE --work /scratch/$USER/mose_tests --one-iteration
python3 test/run_mose_tests.py submit --work /scratch/$USER/mose_tests -p <gpu partition>
python3 test/run_mose_tests.py md5 --work /scratch/$USER/mose_tests > new.md5
python3 test/run_mose_tests.py compare test/reference_1it.md5 new.md5
```

### Checklist

- [ ] The one-iteration table is identical to `test/reference_1it.md5`, or the change explains which entries move and why
- [ ] A change of results comes with a comparison against MOSE on the same input, one step and to the end of a test
- [ ] No new register spills in the production kernels (`ptxinfo` lines of the build log)
- [ ] New subroutines declare `implicit none` and argument intents; no line exceeds 132 characters
- [ ] Documentation updated if user-facing behaviour changed
- [ ] Commit messages use imperative mood ("Add ...", not "Added ...")

## Reporting bugs

Open an issue with the commit, the GPU and the HPC SDK version, the job log, the input file and the smallest case that
reproduces the problem.

## License

MOSE-GPU is released under the [GNU General Public License v3.0](LICENSE). Contributions are accepted under the same
license.

[devguide]: https://open-hydra.github.io/MOSE-GPU/development/contributing/
[testing]: https://open-hydra.github.io/MOSE-GPU/development/testing/
