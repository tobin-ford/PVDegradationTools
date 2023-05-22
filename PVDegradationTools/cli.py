from PVDegradationTools.standards import run_calc_standoff
from PVDegradationTools.humidity import run_calc_rel_humidity
from gaps.cli import CLICommandFromFunction, make_cli

commands = [
    CLICommandFromFunction(
        name="run-standoff", 
        function=run_calc_standoff, 
        split_keys=["project_points"]
    ),
    CLICommandFromFunction(
        name="run-rel_humidity", 
        function=run_calc_rel_humidity, 
        split_keys=["project_points"]
    ),
]
cli = make_cli(commands)


if __name__ == "__main__":
    cli(obj={})