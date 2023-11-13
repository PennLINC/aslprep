#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""ASL preprocessing workflow."""
from aslprep import config

EXITCODE: int = -1


def main():
    """Entry point."""
    import gc
    import sys
    from multiprocessing import Manager, Process
    from os import EX_SOFTWARE
    from pathlib import Path

    from aslprep.cli.parser import parse_args
    from aslprep.cli.workflow import build_workflow
    from aslprep.utils.bids import write_bidsignore, write_derivative_description

    parse_args()

    sentry_sdk = None
    if not config.execution.notrack:
        import sentry_sdk

        from aslprep.utils.sentry import sentry_setup

        sentry_setup()

    # CRITICAL Save the config to a file. This is necessary because the execution graph
    # is built as a separate process to keep the memory footprint low. The most
    # straightforward way to communicate with the child process is via the filesystem.
    config_file = config.execution.work_dir / f"config-{config.execution.run_uuid}.toml"
    config_file.parent.mkdir(exist_ok=True, parents=True)
    config.to_filename(config_file)

    # CRITICAL Call build_workflow(config_file, retval) in a subprocess.
    # Because Python on Linux does not ever free virtual memory (VM), running the
    # workflow construction jailed within a process preempts excessive VM buildup.
    if "pdb" not in config.execution.debug:
        with Manager() as mgr:
            retval = mgr.dict()
            p = Process(target=build_workflow, args=(str(config_file), retval))
            p.start()
            p.join()
            retval = dict(retval.items())  # Convert to base dictionary

            if p.exitcode:
                retval["return_code"] = p.exitcode

    else:
        retval = build_workflow(str(config_file), {})

    global EXITCODE
    EXITCODE = retval.get("return_code", 0)
    aslprep_wf = retval.get("workflow", None)

    # CRITICAL Load the config from the file. This is necessary because the ``build_workflow``
    # function executed constrained in a process may change the config (and thus the global
    # state of ASLPrep).
    config.load(config_file)

    if config.execution.reports_only:
        sys.exit(int(EXITCODE > 0))

    if aslprep_wf and config.execution.write_graph:
        aslprep_wf.write_graph(graph2use="colored", format="svg", simple_form=True)

    EXITCODE = EXITCODE or (aslprep_wf is None) * EX_SOFTWARE
    if EXITCODE != 0:
        sys.exit(EXITCODE)

    # Generate boilerplate
    with Manager() as mgr:
        from aslprep.cli.workflow import build_boilerplate

        p = Process(target=build_boilerplate, args=(str(config_file), aslprep_wf))
        p.start()
        p.join()

    if config.execution.boilerplate_only:
        sys.exit(int(EXITCODE > 0))

    # Clean up master process before running workflow, which may create forks
    gc.collect()

    # Sentry tracking
    if sentry_sdk is not None:
        with sentry_sdk.configure_scope() as scope:
            scope.set_tag("run_uuid", config.execution.run_uuid)
            scope.set_tag("npart", len(config.execution.participant_label))
        sentry_sdk.add_breadcrumb(message="ASLPrep started", level="info")
        sentry_sdk.capture_message("ASLPrep started", level="info")

    config.loggers.workflow.log(
        15,
        "\n".join(["ASLPrep config:"] + [f"\t\t{s}" for s in config.dumps().splitlines()]),
    )
    config.loggers.workflow.log(25, "ASLPrep started!")
    errno = 1  # Default is error exit unless otherwise set
    try:
        aslprep_wf.run(**config.nipype.get_plugin())
    except Exception as e:
        if not config.execution.notrack:
            from aslprep.utils.sentry import process_crashfile

            crashfolders = [
                config.execution.aslprep_dir / f"sub-{s}" / "log" / config.execution.run_uuid
                for s in config.execution.participant_label
            ]
            for crashfolder in crashfolders:
                for crashfile in crashfolder.glob("crash*.*"):
                    process_crashfile(crashfile)

            if "Workflow did not execute cleanly" not in str(e):
                sentry_sdk.capture_exception(e)

        config.loggers.workflow.critical("ASLPrep failed: %s", e)
        raise

    else:
        config.loggers.workflow.log(25, "ASLPrep finished successfully!")
        if not config.execution.notrack:
            success_message = "ASLPrep finished without errors"
            sentry_sdk.add_breadcrumb(message=success_message, level="info")
            sentry_sdk.capture_message(success_message, level="info")

        # Bother users with the boilerplate only iff the workflow went okay.
        boiler_file = config.execution.aslprep_dir / "logs" / "CITATION.md"
        if boiler_file.exists():
            if config.environment.exec_env in (
                "singularity",
                "docker",
                "aslprep-docker",
            ):
                boiler_file = Path("<OUTPUT_PATH>") / boiler_file.relative_to(
                    config.execution.output_dir
                )

            config.loggers.workflow.log(
                25,
                "Works derived from this ASLPrep execution should include the "
                f"boilerplate text found in {boiler_file}.",
            )

        if config.workflow.run_reconall:
            from niworkflows.utils.misc import _copy_any
            from templateflow import api

            dseg_tsv = str(api.get("fsaverage", suffix="dseg", extension=[".tsv"]))
            _copy_any(dseg_tsv, str(config.execution.aslprep_dir / "desc-aseg_dseg.tsv"))
            _copy_any(dseg_tsv, str(config.execution.aslprep_dir / "desc-aparcaseg_dseg.tsv"))

        errno = 0

    finally:
        from fmriprep.reports.core import generate_reports
        from pkg_resources import resource_filename as pkgrf

        # Generate reports phase
        failed_reports = generate_reports(
            config.execution.participant_label,
            config.execution.aslprep_dir,
            config.execution.run_uuid,
            config=pkgrf("aslprep", "data/reports-spec.yml"),
            packagename="aslprep",
        )
        write_derivative_description(
            config.execution.bids_dir,
            config.execution.aslprep_dir,
        )
        write_bidsignore(config.execution.aslprep_dir)

        if failed_reports and not config.execution.notrack:
            sentry_sdk.capture_message(
                f"Report generation failed for {failed_reports} subjects",
                level="error",
            )
        sys.exit(int((errno + failed_reports) > 0))


if __name__ == "__main__":
    raise RuntimeError(
        "aslprep/cli/run.py should not be run directly;\n"
        "Please `pip install` aslprep  and use the `aslprep` command"
    )
