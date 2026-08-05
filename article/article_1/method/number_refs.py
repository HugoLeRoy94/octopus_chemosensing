#!/usr/bin/env python3
"""Bake equation numbers into a LaTeX file for pandoc->docx conversion.

Reads .tex on stdin, writes on stdout:
  - appends  \\qquad(N)  to each numbered display equation (equation/align/...),
    N assigned in document order (matching LaTeX numbering);
  - rewrites \\eqref{eq:X} -> (N) and \\ref{eq:X} -> N;
  - rewrites \\ref{sec:X}/\\ref{tab:X}/\\ref{fig:X} -> the running number, so no
    raw [label] survives (use pandoc --number-sections so headings match).
The source .tex is untouched; this runs on a temp copy inside tex2docx.sh.
"""
import sys, re

text = sys.stdin.read()

# ---- 0. align -> gather (Word/OMML renders the &-alignment '=' as "¿") ----
# Each alignment '&=' becomes a plain '=', each column-separator '&' becomes a
# line break, so a 2-row/2-column align becomes stacked, centered lines that
# render correctly. The source .tex keeps its align (this runs on a temp copy).
def align_to_gather(m):
    body = m.group(1).replace('&=', '=').replace('&', r'\\')
    return r'\begin{gather}' + body + r'\end{gather}'

text = re.sub(r'\\begin\{align\*?\}(.*?)\\end\{align\*?\}',
              align_to_gather, text, flags=re.DOTALL)

# ---- 1. number display equations, insert \qquad(N) before each \label{eq:} ----
eqnum, counter = {}, [0]

def do_env(m):
    env, body = m.group(1), m.group(2)
    if env in ('align', 'gather', 'eqnarray', 'multline', 'flalign'):
        segs = body.split('\\\\')
        for k, seg in enumerate(segs):
            labs = re.findall(r'\\label\{(eq:[^}]*)\}', seg)
            if labs:
                counter[0] += 1
                for lb in labs:
                    eqnum[lb] = counter[0]
                segs[k] = re.sub(r'(\\label\{eq:)',
                                 r'\\qquad(%d)\1' % counter[0], seg, count=1)
        body = '\\\\'.join(segs)
    else:  # equation
        labs = re.findall(r'\\label\{(eq:[^}]*)\}', body)
        if labs:
            counter[0] += 1
            for lb in labs:
                eqnum[lb] = counter[0]
            body = re.sub(r'(\s*\\label\{eq:)',
                          r' \\qquad(%d)\1' % counter[0], body, count=1)
    return '\\begin{%s}%s\\end{%s}' % (env, body, env)

text = re.sub(r'\\begin\{(equation|align|gather|multline|eqnarray|flalign)\}(.*?)\\end\{\1\}',
              do_env, text, flags=re.DOTALL)

# ---- 2. map sec/tab/fig labels to running numbers (document order) ----
secnum = {lb: i for i, lb in enumerate(re.findall(r'\\label\{(sec:[^}]*)\}', text), 1)}
tabnum = {lb: i for i, lb in enumerate(re.findall(r'\\label\{(tab:[^}]*)\}', text), 1)}
fignum = {lb: i for i, lb in enumerate(re.findall(r'\\label\{(fig:[^}]*)\}', text), 1)}

# ---- 3. rewrite all internal references ----
def repl(m):
    cmd, lb = m.group(1), m.group(2)
    if lb in eqnum:
        return ('(%d)' % eqnum[lb]) if cmd == 'eqref' else str(eqnum[lb])
    for d in (secnum, tabnum, fignum):
        if lb in d:
            return str(d[lb])
    return m.group(0)

text = re.sub(r'\\(eqref|ref)\{([^}]*)\}', repl, text)

sys.stdout.write(text)
