#!/usr/bin/env python3
"""The regression table of MOSE-GPU on the test cases of MOSE.

The test cases of MOSE that run on the GPU are kept in this folder, ready to run, at their path in test/ of MOSE
(1D/Sod79, 2D/euler/chang-nozzle, ...). Stage them, submit them with the `mose` launcher, collect the md5 of the final
fields, and compare two tables.

  run_mose_tests.py stage   --work WORK [--one-iteration] [--mose MOSE --swbli-seed DIR] [TEST ...]
  run_mose_tests.py submit  --work WORK [-p PARTITION] [-t TIME] [--mem MEM] [TEST ...]
  run_mose_tests.py md5     --work WORK [TEST ...]
  run_mose_tests.py compare REFERENCE NEW
  run_mose_tests.py import  --mose MOSE [TEST ...]

TEST is a path under test/, e.g. 1D/Sod79, or its table name 1D_Sod79; by default the 23 cases of
test/reference_1it.md5.

stage copies the cases into WORK, one folder per table name. --one-iteration sets iter-threshold = 1: the regression
table. swbli SA and SST are not kept here: MOSE-GPU has no grid sequencing, so the case restarts on the fine grid
(levels = 1, newrun = false) from the field of the coarse levels computed by MOSE. They are staged from a MOSE checkout,
--mose MOSE --swbli-seed DIR with DIR/SA/field.tec and DIR/SST/field.tec; SST also gets point-implicit = true. Without
these two options the two cases are skipped.

import writes the cases of a MOSE checkout into this folder (the copies here come from MOSE b29130b), the way the
validation against MOSE staged them: the files of INPUT/ and MESH/ of the MOSE test, then those the `solve` block of its
MOSE.sh links, as real files; input.ini without carriage returns. As in the validation:
  - Ein91, Noh87, Tor99: the 400-cell grid (ic_x4, bc_x4), like Sod79;
  - isentropic-vortex, mms: the 64 x 64 grid (ic-64, bc-64); mms runs without its manufactured source;
  - gresho-vortex: the solve argument 0.01 (Mach number);
  - an iteration cap: 5000 for the nozzles, 2000 for pressure-centrifugal.
The initial fields and boundary files of the grids a case does not use (ic_x8.tec, bc-128.txt, ...) are left out.

Python 3.6 or newer, standard library only.
"""
import argparse
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
MOSE_LAUNCHER = os.path.join(os.path.dirname(HERE), 'mose')
REFERENCE = os.path.join(HERE, 'reference_1it.md5')
GRID_X4 = ('1D/Ein91', '1D/Noh87', '1D/Tor99')
GRID_64 = ('2D/euler/isentropic-vortex', '2D/viscous/mms')
CAPS = {'2D/euler/chang-nozzle': 5000, '2D/euler/rocket-nozzle/frozen': 5000, '3D/pressure-centrifugal': 2000}
SOLVE_ARG = {'2D/euler/gresho-vortex': '0.01'}
SWBLI = ('2D/viscous/swbli/SA', '2D/viscous/swbli/SST')
GRID_VARIANT = re.compile(r'^(ic|bc)[-_].+\.(tec|txt)$')      # ic_x8.tec, bc-128.txt, ic-coarse.tec, ...
RUN_FILES = shutil.ignore_patterns('OUTPUT', 'mose_*.out', 'mose_run.sbatch')


def name_of(test):
    return test.replace('/', '_')


# the cases of the regression table, in the order of test/reference_1it.md5
TESTS = ('1D/diff1', '1D/diff2', '1D/Ein91', '1D/Noh87', '1D/Shu-Osher', '1D/Sod79', '1D/Tor99',
         '2D/euler/chang-nozzle', '2D/euler/gresho-vortex', '2D/euler/hypersonic-cylinder', '2D/euler/isentropic-vortex',
         '2D/euler/oblique-shock', '2D/euler/prandtl-meyer', '2D/euler/ramp-channel', '2D/euler/rocket-nozzle/frozen',
         '2D/euler/supersonic-forward-step', '2D/viscous/flat-plate-laminar', '2D/viscous/flat-plate-turbulent/SA',
         '2D/viscous/flat-plate-turbulent/SST', '2D/viscous/mms', '2D/viscous/swbli/SA', '2D/viscous/swbli/SST',
         '3D/pressure-centrifugal')


def test_path(t):
    """A case given as a path (1D/Sod79) or as its table name (1D_Sod79)."""
    if '/' in t:
        return t
    for p in TESTS:
        if name_of(p) == t:
            return p
    raise SystemExit('%s: unknown case name; give the path under test/ of MOSE' % t)


def read_table(path):
    table = {}
    for line in open(path):
        w = line.split()
        if len(w) >= 2 and not line.startswith('#'):
            table[w[0]] = w[1]
    return table


def set_key(ini, section, key, value):
    """Replace key in [section], or add it at the top of the section; add the section when absent."""
    lines = open(ini).read().split('\n')
    out, cur, done, head = [], None, False, None
    for line in lines:
        m = re.match(r'\s*\[(.+)\]', line)
        if m:
            cur = m.group(1).strip()
            if cur == section:
                head = len(out)
        elif cur == section and re.match(r'\s*' + re.escape(key) + r'\s*=', line):
            line, done = '%s = %s' % (key, value), True
        out.append(line)
    if not done:
        if head is None:
            out = ['[%s]' % section, '%s = %s' % (key, value), ''] + out
        else:
            out.insert(head + 1, '%s = %s' % (key, value))
    open(ini, 'w').write('\n'.join(out))


def get_key(ini, section, key):
    cur = None
    for line in open(ini):
        m = re.match(r'\s*\[(.+)\]', line)
        if m:
            cur = m.group(1).strip()
            continue
        m = re.match(r'\s*' + re.escape(key) + r'\s*=\s*(.*)', line)
        if cur == section and m:
            return m.group(1).split(';')[0].strip()
    return ''


def link(target, path):
    if os.path.lexists(path):
        os.remove(path)
    os.symlink(target, path)


def solve_links(src, dst, arg):
    """The links the `solve` block of the test's MOSE.sh makes before it runs, made in the staged copy."""
    script = os.path.join(src, 'MOSE.sh')
    if not os.path.isfile(script):
        return
    lines = open(script).read().splitlines()
    start = next((i for i, l in enumerate(lines) if re.search(r'\$1 == solve', l)), None)
    if start is None:
        return
    rel = '.'
    for l in lines[start + 1:]:
        s = l.split('#')[0].strip()
        if s.startswith('mkdir -p OUTPUT'):
            break
        w = s.replace('$2', arg).split()
        if not w:
            continue
        if w[0] == 'cd':
            rel = os.path.normpath(os.path.join(rel, w[1]))
        elif w[0] == 'ln':
            a, b = [x for x in w[1:] if not x.startswith('-')]
            target = os.path.realpath(os.path.join(src, rel, a))
            if not os.path.isfile(target):
                raise SystemExit('%s: MOSE.sh solve links %s, which does not exist' % (src, a))
            link(target, os.path.join(dst, rel, b))


def from_mose(t, src, dst, swbli_seed=None):
    """Test t of MOSE (folder src) staged into dst as the validation staged it: links to the files of MOSE."""
    if os.path.isdir(dst):
        shutil.rmtree(dst)
    for d in ('INPUT', 'MESH', 'OUTPUT'):
        os.makedirs(os.path.join(dst, d))
    for d in ('INPUT', 'MESH'):
        sd = os.path.join(src, d)
        if os.path.isdir(sd):
            for f in sorted(os.listdir(sd)):
                p = os.path.join(sd, f)
                if os.path.isfile(p):
                    link(os.path.realpath(p), os.path.join(dst, d, f))
    solve_links(src, dst, SOLVE_ARG.get(t, ''))
    if t in GRID_X4 or t in GRID_64:
        tag = '_x4' if t in GRID_X4 else '-64'
        link(os.path.realpath(os.path.join(src, 'INPUT', 'ic%s.tec' % tag)), os.path.join(dst, 'INPUT', 'ic.tec'))
        link(os.path.realpath(os.path.join(src, 'INPUT', 'bc%s.txt' % tag)), os.path.join(dst, 'INPUT', 'bc.txt'))
    ini = os.path.join(dst, 'input.ini')
    open(ini, 'w').write(open(os.path.join(src, 'input.ini')).read().replace('\r', ''))
    if t in CAPS:
        old = get_key(ini, 'MOSE-Parameters', 'iter-threshold')
        if not old or int(old) > CAPS[t]:
            set_key(ini, 'MOSE-Parameters', 'iter-threshold', str(CAPS[t]))
    if t in SWBLI:
        seed = os.path.join(os.path.abspath(swbli_seed), os.path.basename(t), 'field.tec')
        if not os.path.isfile(seed):
            raise SystemExit('%s: no seed field %s' % (t, seed))
        shutil.copy(seed, os.path.join(dst, 'OUTPUT', 'field.tec'))
        set_key(ini, 'MOSE-Multigrid', 'levels', '1')
        set_key(ini, 'MOSE-Parameters', 'newrun', 'false')
        set_key(ini, 'MOSE-Parameters', 'iter-threshold', '20000')
        if t.endswith('SST'):
            set_key(ini, 'MOSE-Turbulence', 'point-implicit', 'true')


def stage(args):
    tests = [test_path(t) for t in args.tests] if args.tests else list(TESTS)
    for t in tests:
        src, dst = os.path.join(HERE, t), os.path.join(os.path.abspath(args.work), name_of(t))
        if os.path.isfile(os.path.join(src, 'input.ini')):
            if os.path.isdir(dst):
                shutil.rmtree(dst)
            shutil.copytree(src, dst, ignore=RUN_FILES)
            for d in ('INPUT', 'MESH', 'OUTPUT'):
                os.makedirs(os.path.join(dst, d), exist_ok=True)
        elif t in SWBLI:
            if not (args.mose and args.swbli_seed):
                print('%-45s skipped: needs --mose and --swbli-seed' % t)
                continue
            src = os.path.join(os.path.abspath(args.mose), 'test', t)
            if not os.path.isdir(src):
                print('%-45s no such test in %s' % (t, os.path.dirname(src)))
                continue
            from_mose(t, src, dst, args.swbli_seed)
        else:
            print('%-45s no such case in %s' % (t, HERE))
            continue
        if args.one_iteration:
            set_key(os.path.join(dst, 'input.ini'), 'MOSE-Parameters', 'iter-threshold', '1')
        print('%-45s staged' % t)


def import_cases(args):
    """The cases of a MOSE checkout written into this folder, with real files in place of the links of from_mose."""
    tests = [test_path(t) for t in args.tests] if args.tests else [t for t in TESTS if t not in SWBLI]
    root = os.path.join(os.path.abspath(args.mose), 'test')
    tmp = tempfile.mkdtemp()
    try:
        for t in tests:
            src = os.path.join(root, t)
            if t in SWBLI:
                print('%-45s not imported: staged from MOSE with --mose and --swbli-seed' % t)
                continue
            if not os.path.isdir(src):
                print('%-45s no such test in %s' % (t, root))
                continue
            work = os.path.join(tmp, name_of(t))
            from_mose(t, src, work)
            dst = os.path.join(HERE, t)
            if os.path.isdir(dst):
                shutil.rmtree(dst)
            os.makedirs(dst)
            ini = open(os.path.join(work, 'input.ini')).read()
            shutil.copyfile(os.path.join(work, 'input.ini'), os.path.join(dst, 'input.ini'))
            left = []
            for d in ('INPUT', 'MESH'):
                names = sorted(os.listdir(os.path.join(work, d)))
                for f in names:
                    m = GRID_VARIANT.match(f)
                    if d == 'INPUT' and m and '%s.%s' % m.groups() in names and f not in ini:
                        left.append(f)      # the field or boundary file of a grid the case does not use
                        continue
                    os.makedirs(os.path.join(dst, d), exist_ok=True)
                    shutil.copyfile(os.path.join(work, d, f), os.path.join(dst, d, f))
            print('%-45s imported%s' % (t, ', left out ' + ' '.join(left) if left else ''))
    finally:
        shutil.rmtree(tmp)


def submit(args):
    work = os.path.abspath(args.work)
    names = [name_of(t) if '/' in t else t for t in args.tests] if args.tests else sorted(os.listdir(work))
    for n in names:
        d = os.path.join(work, n)
        if not os.path.isfile(os.path.join(d, 'input.ini')):
            continue
        cmd = [sys.executable, MOSE_LAUNCHER, 'submit', '-n', '1', '-t', args.time, '-J', 't1_' + n, '.']
        if args.partition:
            cmd[3:3] = ['-p', args.partition]
        if args.mem:
            cmd[3:3] = ['--mem', args.mem]
        p = subprocess.run(cmd, cwd=d, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        out = p.stdout.decode(errors='replace').splitlines()
        lines = [l for l in out if '->' in l or l.startswith('ERROR')] or out[-1:]
        for l in lines:
            print('%-45s %s' % (n, l.strip()[:120]))


def md5(args):
    work = os.path.abspath(args.work)
    names = [name_of(t) if '/' in t else t for t in args.tests] if args.tests else sorted(os.listdir(work))
    for n in names:
        d = os.path.join(work, n)
        if not os.path.isdir(d):
            continue
        logs = sorted(f for f in os.listdir(d) if f.startswith('mose_') and f.endswith('.out'))
        text = open(os.path.join(d, logs[-1]), errors='replace').read() if logs else ''
        if re.search(r'^\[mose\] end .*rc=0', text, re.M):
            h = hashlib.md5(open(os.path.join(d, 'OUTPUT', 'field.tec'), 'rb').read()).hexdigest()[:8]
        elif re.search(r'^\[mose\] end', text, re.M):
            h = 'FAILED'
        else:
            h = 'RUNNING'
        print('%-45s %s' % (n, h))


def compare(args):
    ref, new = read_table(args.reference), read_table(args.new)
    same = diff = pending = 0
    for n in sorted(set(ref) | set(new)):
        a, b = ref.get(n), new.get(n)
        if a is None or b is None:
            print('%-45s %s' % (n, 'not in the reference' if a is None else 'not run'))
        elif b in ('RUNNING', 'FAILED'):
            pending += 1
            print('%-45s %s' % (n, 'not finished yet' if b == 'RUNNING' else 'FAILED: read its mose_*.out'))
        elif a == b:
            same += 1
        else:
            diff += 1
            print('%-45s DIFFERENT: %s -> %s' % (n, a, b))
    print('identical: %d   different: %d   not finished or failed: %d' % (same, diff, pending))
    return 1 if diff or pending else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    sub = ap.add_subparsers(dest='cmd')
    p = sub.add_parser('stage')
    p.add_argument('--work', required=True, help='where to stage the cases')
    p.add_argument('--one-iteration', action='store_true', help='iter-threshold = 1: the regression table')
    p.add_argument('--mose', help='a checkout of MOSE (open-hydra/MOSE), for the two swbli cases')
    p.add_argument('--swbli-seed', help='DIR/SA/field.tec and DIR/SST/field.tec, the coarse-level fields of MOSE')
    p.add_argument('tests', nargs='*')
    p = sub.add_parser('submit')
    p.add_argument('--work', required=True)
    p.add_argument('-p', '--partition', default='')
    p.add_argument('-t', '--time', default='00:20:00')
    p.add_argument('--mem', default='', help='memory per job, passed to mose submit')
    p.add_argument('tests', nargs='*')
    p = sub.add_parser('md5')
    p.add_argument('--work', required=True)
    p.add_argument('tests', nargs='*')
    p = sub.add_parser('compare')
    p.add_argument('reference')
    p.add_argument('new')
    p = sub.add_parser('import')
    p.add_argument('--mose', required=True, help='a checkout of MOSE (open-hydra/MOSE)')
    p.add_argument('tests', nargs='*')
    args = ap.parse_args()
    if args.cmd == 'stage':
        stage(args)
    elif args.cmd == 'submit':
        submit(args)
    elif args.cmd == 'md5':
        md5(args)
    elif args.cmd == 'compare':
        sys.exit(compare(args))
    elif args.cmd == 'import':
        import_cases(args)
    else:
        ap.print_help()


if __name__ == '__main__':
    main()
