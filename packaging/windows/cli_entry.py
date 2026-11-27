import sys

from hischeck.__main__ import main


if __name__ == "__main__":
    # The desktop runner reads both pipes as UTF-8, including paths and errors.
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    raise SystemExit(main())
