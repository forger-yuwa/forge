"""cwd が指定の run の forge を止める (SIGTERM → 180 s 待って SIGKILL → 終了を確かめる)。止まれば終了コード 0。usage: python3 m9_stop.py <run>"""
import os, signal, subprocess, sys, time
run = os.path.realpath(sys.argv[1])
def pids():
    r = []
    for p in subprocess.run(["pgrep", "-x", "forge"], capture_output=True, text=True).stdout.split():
        try:
            if os.path.realpath(f"/proc/{p}/cwd") == run: r.append(int(p))
        except OSError: pass
    return r
for sig, wait in ((signal.SIGTERM, 180), (signal.SIGKILL, 60)):
    for p in pids():
        try: os.kill(p, sig)
        except OSError: pass
    t0 = time.time()
    while pids() and time.time() - t0 < wait: time.sleep(3)
sys.exit(0 if not pids() else 1)
