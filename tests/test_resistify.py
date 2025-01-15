import subprocess
import filecmp
import os


def check_results(expected_directory, actual_directory):
    expected_results_file = os.path.join(expected_directory, "results.tsv")
    actual_results_file = os.path.join(actual_directory, "results.tsv")

    assert filecmp.cmp(expected_results_file, actual_results_file), (
        "Results file mismatch"
    )


# def test_resistify_download_models():
#    # Define the models directory
#    models_dir = "models"
#
#    # Run the download_models command
#    result = subprocess.run(
#        ["resistify", "download_models", models_dir],
#        capture_output=True,
#        text=True,
#    )
#
#    # Check the exit code
#    assert result.returncode == 0, f"Download models failed: {result.stderr}"
#
#    # Check that the models directory exists
#    assert os.path.isdir(models_dir), "Models directory not created"


def test_resistify_nlr():
    # Define paths to input and output files
    input_file = "tests/data/zar1.fa"
    expected_directory = "tests/data/nlr_expected"
    actual_directory = "tests/output/nlr_actual"

    # Run the Resistify command for NLR
    result = subprocess.run(
        ["resistify", "nlr", input_file, "-o", actual_directory],
        capture_output=True,
        text=True,
    )

    # Check the exit code
    assert result.returncode == 0, f"Command failed: {result.stderr}"

    check_results(expected_directory, actual_directory)


def test_resistify_nlr_retain():
    # Define paths to input and output files
    input_file = "tests/data/zar1.fa"
    expected_directory = "tests/data/nlr_retain_expected"
    actual_directory = "tests/output/nlr_retain_actual"

    # Run the Resistify command for NLR
    result = subprocess.run(
        ["resistify", "nlr", "--retain", input_file, "-o", actual_directory],
        capture_output=True,
        text=True,
    )

    # Check the exit code
    assert result.returncode == 0, f"Command failed: {result.stderr}"

    check_results(expected_directory, actual_directory)


def test_resistify_nlr_coconat():
    # Define paths to input and output files
    input_file = "tests/data/zar1.fa"
    expected_output = "tests/data/nlr_coconat_expected"
    actual_output = "tests/output/nlr_coconat_actual"

    # Run the Resistify command for NLR
    result = subprocess.run(
        ["resistify", "nlr", "--coconat", input_file, "-o", actual_output],
        capture_output=True,
        text=True,
    )

    # Check the exit code
    assert result.returncode == 0, f"Command failed: {result.stderr}"

    check_results(expected_output, actual_output)


def test_resistify_nlr_retain_coconat():
    # Define paths to input and output files
    input_file = "tests/data/zar1.fa"
    expected_directory = "tests/data/nlr_retain_coconat_expected"
    actual_directory = "tests/output/nlr_retain_coconat_actual"

    # Run the Resistify command for NLR
    result = subprocess.run(
        [
            "resistify",
            "nlr",
            "--retain",
            "--coconat",
            input_file,
            "-o",
            actual_directory,
        ],
        capture_output=True,
        text=True,
    )

    # Check the exit code
    assert result.returncode == 0, f"Command failed: {result.stderr}"

    check_results(expected_directory, actual_directory)


def test_resistify_prr():
    # Define paths to input and output files
    input_file = "tests/data/fls2.fa"
    expected_output = "tests/data/prr_expected"
    actual_output = "tests/output/prr_actual"

    # Run the Resistify PRR command
    result = subprocess.run(
        ["resistify", "prr", input_file, "-o", actual_output],
        capture_output=True,
        text=True,
    )

    # Check the exit code
    assert result.returncode == 0, f"PRR command failed: {result.stderr}"

    check_results(expected_output, actual_output)
