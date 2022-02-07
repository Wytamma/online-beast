import os
from typer.testing import CliRunner

from online_beast.main import app

runner = CliRunner()


def test_app():
    result = runner.invoke(
        app,
        [
            "data/testGTR.xml",
            "data/samples.fasta",
            "--state-file",
            "data/testGTR.xml.state",
            "--output",
            "testGTR.xml",
        ],
    )
    assert result.exit_code == 0
    os.remove("testGTR.xml")
    os.remove("testGTR.xml.state")
