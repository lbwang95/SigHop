#!/usr/bin/env python3
"""Forward job output and terminate its process group on a mismatch or signal."""
import os, re, signal, subprocess, sys

child=None

def stop_child():
    if child is None or child.poll() is not None:return
    try:os.killpg(child.pid,signal.SIGTERM)
    except ProcessLookupError:return
    try:child.wait(timeout=3)
    except subprocess.TimeoutExpired:
        try:os.killpg(child.pid,signal.SIGKILL)
        except ProcessLookupError:pass
        child.wait()

def signal_stop(sig,frame):
    stop_child()
    raise SystemExit(128+sig)

def main():
    global child
    signal.signal(signal.SIGTERM,signal_stop);signal.signal(signal.SIGINT,signal_stop)
    bad=re.compile(r'\bmismatches\s*=\s*([1-9][0-9]*)')
    oom=False;mismatch=False
    try:
        child=subprocess.Popen(['bash','-c',sys.argv[1]],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,
                               text=True,encoding='utf-8',errors='replace',start_new_session=True)
        for line in child.stdout:
            sys.stdout.write(line);sys.stdout.flush()
            oom=oom or any(s in line.lower() for s in ('std::bad_alloc','out of memory','cannot allocate memory','memoryerror'))
            if bad.search(line):
                mismatch=True;stop_child();break
        code=child.wait()
        print('# child_exit_code=%d' % code,flush=True)
        if mismatch:
            print('ERROR: nonzero mismatches; stopped child group',flush=True);return 1
        if code and oom:
            print('# explicit allocation failure classified as MO',flush=True);return 4
        return code if code>=0 else 128-code
    finally:stop_child()

if __name__=='__main__':sys.exit(main())
