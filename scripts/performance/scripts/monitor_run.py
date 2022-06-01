import platform
import subprocess

import typer


def main(script_path: str, tree_path: str, out_path: str):
    exe = 'gtime' if platform.system() == 'Darwin' else '/usr/bin/time'
    args = [exe, '--verbose', 'python', script_path, tree_path]
    proc = subprocess.Popen(args, stderr=subprocess.PIPE, encoding='utf-8')
    stdout, stderr = proc.communicate()
    if proc.returncode == 0:
        with open(out_path, 'w') as fh:
            fh.write(stderr)
    else:
        raise Exception(f'Non-zero return code: {stderr}')
    return


if __name__ == "__main__":
    typer.run(main)
