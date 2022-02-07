import os
from online_beast import __version__
from typer.testing import CliRunner

from online_beast.main import app

runner = CliRunner()


def test_version():
    assert __version__ == "0.1.0"


def test_app():
    result = runner.invoke(
        app,
        [
            "data/testGTR.xml",
            "data/samples.fasta",
            "--state-file",
            "data/testGTR.xml.state",
            "--xml-output",
            "testGTR.xml",
            "--sate-output",
            "testGTR.xml.state",
        ],
    )
    assert result.exit_code == 0
    os.remove("testGTR.xml")
    os.remove("testGTR.xml.state")
