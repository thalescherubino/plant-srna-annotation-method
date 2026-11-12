# Step 4: Environments, compute and the run queue

[Back to the skill map](../SKILL.md) · Previous: [03_library_processing_qc.md](03_library_processing_qc.md) · Next: [05_mirna_prediction.md](05_mirna_prediction.md)

## Environments

Each tool gets its own conda environment, defined in the repository's `envs/` folder; step-by-step installation is in `docs/setup.md`.

| Environment | Tools |
|---|---|
| `plant-srna-phasis` | Phasis 2.8.x-2.9.x (installed from its source checkout with the pip constraints file; record the exact version and commit), HISAT2, samtools |
| `plant-srna-shortstack` | ShortStack 4.1.2 |
| `plant-srna-mirador` | miRador (Python 3.12 port), Bowtie 1, BLAST, EMBOSS, ViennaRNA |
| `plant-srna-mirscore` | miRScore 0.3.7, Bowtie 1, samtools, ViennaRNA |
| `plant-srna-structural` | tRNAscan-SE, Infernal, Bowtie 1 |
| `plant-srna-reads` | sra-tools, cutadapt, Bowtie 1, pigz (SRA route) |

On a shared server:
- Install a **private** Miniforge in your own work directory, with its own `.condarc`. Never run `conda init`, never install into a shared conda.
- Run tools without activating: `mamba run -n <env> <command>`.
- Record versions and source commits (`phasis -version`, `ShortStack --version`, `git rev-parse HEAD`) in a `software_versions.txt` next to the results.
- Phasis pins NumPy and scikit-learn but not pandas or matplotlib; without the constraints file, pip pulls the newest pandas (3.x), which was not tested.

## Choosing thread counts

- Check the load first (`uptime`, `top`). `nproc` can report 1 when `OMP_NUM_THREADS=1`; use `getconf _NPROCESSORS_ONLN` for the real CPU count.
- Guideline: about 32 cores per 100 libraries for heavy steps, one heavy job at a time on a shared server. The 24-nt Phasis runs benefit from more cores (we used 64 when the server was under half load). These are guidelines, not requirements, and Phasis is being optimized, so expect them to change.
- Phasis treats `-cores` as a ceiling: heavy stages use up to 70% of it, and many stages use far fewer.
- Ask before running two heavy steps in parallel.

## Run long jobs detached, through a queue

Most steps take hours; the whole species takes 2-4 days. Never tie a run to a laptop or an editor session.

- Start the queue in `tmux` (or `setsid nohup`), with output going to a log file.
- Each step appends `start <name> <time>`, then `done` or `failed rc=<code>` to a status file (`runs/QUEUE_STATUS`). A step marked `done` is skipped when the queue is re-run, so a failed queue can simply be restarted after a fix.
- Each step writes its own log (`run.<step>.log`) in its output folder.
- Make every step **fail loudly**. A shell function's exit code is that of its last command. If an earlier `cp` fails silently, the step is logged `done` with missing output; we lost a miRScore run this way. End functions with an explicit check or `|| return 1`.

Queue order used for tomato (one 32-core job at a time):

| # | Step | Tomato time |
|---|---|---|
| 1 | exact-mapping QC | 0.3 h |
| 2 | structural RNA scan | 6.8 h |
| 3 | Phasis 21 nt, individual libraries (builds the shared index) | 4.3 h |
| 4 | ShortStack miRNA discovery | 2.6 h |
| 5 | miRador | about 1 h (rice) |
| 6 | Phasis 24 nt, individual | more than 28 h |
| 7 | ShortStack with known miRNAs | 3.0 h |
| 8 | expand libraries to uncompressed FASTA | 4 min |
| 9 | miRScore inputs (+ hairpinHelper for miRBase) | 1 min |
| 10 | miRScore × 3 sets (8 chunks × 4 threads) | 3.8-5 h each |
| 11 | Phasis 21 and 24 nt, pooled | rice: pooled 24 took 28.5 h |
| 12 | `phasis-compare` per phase | minutes |

Rice (374 Mb genome) took about 60 hours at 32 cores in total. Time grows with genome size and is dominated by the 24-nt Phasis runs.

## Changing a queue while it runs

- **Never edit a running bash script in place.** Bash reads the file as it goes, so an edit shifts what it reads next. Write a new script instead.
- To run a later step early (for example the miRNA steps while Phasis is busy), pause the queue controller with `kill -STOP <pid>`, run the steps from a side script that uses the same status file, then `kill -CONT <pid>`. Children of the paused controller (Phasis itself) keep running. Make the side script resume the controller on exit, even on failure.
- `pkill -f` and `pgrep -f` patterns can match your own shell. Kill by PID.
- Ask before deleting large intermediates (expanded FASTA, Phasis `processed_libraries/`).

## Monitoring

- Phasis logs progress bars with elapsed time per stage; a small sampler logging CPU and memory every minute helps locate slow stages.
- A log that stops changing is not proof of a hang. Check whether a child process is using CPU (`ps`, `pstree -ap <pid>`). miRScore's rescue step, for example, runs one `samtools view` per library per hairpin and prints only in bursts.
