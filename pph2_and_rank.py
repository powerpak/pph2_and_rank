#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""\
pph2_and_rank.py

Submits a list of genes with SNPs to PolyPhen2 and orders them by
a score measuring the overall deleteriousness of the mutations, 
normalized for protein length.

LICENSE

Written by Ted Pak for the Roth Laboratory (http://llama.mshri.on.ca).

Released under an MIT license; please see MIT-LICENSE.txt.
"""

USAGE = """\
Usage: %s [-h|--help] [-q|--quiet] [-s|--sid SID] [-o|--output OUTPUT]
       MUT_LIST_1 [MUT_LIST_2 ...]
Options:
  -h|--help             Displays this message.
  -o|--output OUTPUT    Prints to file named OUTPUT instead of standard output.
  -q|--quiet            Do not print status messages to standard error.
  -s|--sid SID          Specify a pre-existing GGI SID.
                        Use this to continue from a prior run.
See README.mdown for complete details."""

import sys, urllib, urllib2, csv, time, math, os
from BeautifulSoup import BeautifulSoup
import re
import sqlite3
import socket
socket.setdefaulttimeout(10)

pph2_url = 'http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi2.cgi'
pph2_result_url = 'http://genetics.bwh.harvard.edu/ggi/pph2/%s/1/pph2-full.txt'
pph2_track_url = 'http://genetics.bwh.harvard.edu/cgi-bin/ggi/ggi2.cgi?_ggi_project=PPHWeb2'\
    + '&sid=%s&_ggi_target_manage=Refresh&_ggi_origin=manage'

# If the list is too big to fit in memory, change this to a filename.
db = sqlite3.connect(':memory:')
c = db.cursor()
d = db.cursor()

accid_cache = {}
quiet = False
steps = ['(1/7) Validating input', '(2/7) Mapping genomic SNPs', 
    '(3/7) Collecting output', '(4/7) Building MSA and annotating proteins',
    '(5/7) Collecting output', '(6/7) Predicting',
    '(7/7) Generating reports']

# Attempts to get the terminal width; will only succeed on a *NIX
try:
    term_columns = map(int, os.popen('stty size', 'r').read().split())[1]
except Exception:
    term_columns = 80

def setup_db():
    """Create the SQLite3 database that will hold interim results."""
    c.execute('''CREATE TABLE t (
      id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
      gene TEXT DEFAULT NULL,
      accid TEXT DEFAULT NULL,
      seqlen INTEGER DEFAULT NULL,
      seqpos INTEGER DEFAULT NULL,
      aa1 TEXT DEFAULT NULL,
      aa2 TEXT DEFAULT NULL,
      pph2_score REAL
    )''')
    c.execute('CREATE INDEX idx_gene ON t (gene ASC)')
    c.execute('CREATE INDEX idx_accid ON t (accid ASC)')

def write_status(msg, num=None, denom=None):
    """ Writes a progress message to standard error.
        num and denom can be used to add a progress bar."""
    if quiet: return
    if msg is True:
        return sys.stderr.write(' => Done.\n')
    if num is not None: 
        sys.stderr.write('\r' + ''.join([' '] * term_columns) + '\r')
        if denom is not None: msg  = "%s => %s" % (msg, progress(num, denom))
    sys.stderr.write(msg)

def progress(curr, finish, width=10):
    """Creates a progress bar using unicode characters; the terminal must support unicode."""
    ticks = ['▂','▃','▄','▅','▆','▇']
    finish = max(finish, 1)
    pct = max(curr, 0) / float(finish)
    bartip, barlen = math.modf(pct * width)
    bartip = ticks[int(len(ticks) * bartip)] if pct < 1.0 else ''
    bar = "%s%s%s" % (ticks[-1] * int(barlen), bartip, ' ' * (width - int(barlen) - 1))
    prog = "%s%% [ %s ] %d/%d" % (str(int(pct*100)).rjust(3), bar, curr, finish)
    return prog

def spin(seconds):
    """Count down a number of seconds, printing the countdown to standard error."""
    for x in range(seconds):
        padded_secs = str(seconds - x).rjust(len(str(seconds)))
        msg = ' (%s)' % padded_secs
        sys.stderr.write(msg)
        time.sleep(1)
        sys.stderr.write('\b' * len(msg))
    
def get_seqlen_and_gene(accid):
    """Gets the sequence length and gene name for a SwissProt/UniProt accession ID."""
    cached = accid_cache.get(accid)
    if cached is not None:
        return cached
    result = [None, None]
    try:
        f = urllib.urlopen("http://www.uniprot.org/uniprot/%s.txt" % accid)
        lines = f.readlines()
    except (socket.timeout, IOError):
        raise RemoteException("UniProt didn't respond to an HTTP request.  That's odd.")
    for line in lines:
        match = re.match(r'^SQ\s+SEQUENCE\s+(\d+)\s*AA;', line)
        if match: result[0] = match.group(1)
        match = re.match(r'^GN   Name=([^;]+);', line)
        if match: result[1] = match.group(1)
        if result[0] is not None and result[1] is not None:
            accid_cache[accid] = result
            break
    return result

def submit_to_polyphen2(batch, humvar=True, hg19=False):
    """ Submits a batch file of accids + SNPs to PolyPhen2 and returns 
        the session ID for the job."""
    params = urllib.urlencode({
        '_ggi_project': 'PPHWeb2', 
        '_ggi_origin': 'query', 
        '_ggi_target_pipeline': '1',
        'MODELNAME': 'HumVar' if humvar else 'HumDiv',
        'UCSCDB': 'hg19' if hg19 else 'hg18',
        'SNPFUNC': 'm',
        'NOTIFYME': '',
        'SNPFILTER': '0',
        '_ggi_batch': batch
    })
    doc = None
    while doc is None:
        try:
            response = urllib2.urlopen(pph2_url, params)
            doc = response.read()
        except (socket.timeout, IOError): pass
    soup = BeautifulSoup(doc)
    sid_input = soup.find('input', {'name': 'sid'})
    if sid_input is None:
        print doc
        raise RemoteException("GGI returned a weird page without a SID.")
    sid = sid_input['value']
    return sid

def poll_for_polyphen2_results(sid):
    """ Polls PolyPhen2's GGI web interface for updates on the progress of the job.
        Once the job has completed the full result file is returned. """
    curr_step = -1
    max_tries = 10
    tries = 0
    wait_msg = "Waiting for PolyPhen2 results => %s"
    done_msg = " => Done.\n"
    while True:
        params = urllib.urlencode({
            '_ggi_project': 'PPHWeb2', 
            '_ggi_origin': 'manage', 
            '_ggi_target_manage': 'Refresh',
            'sid': sid
        })
        doc = None
        while doc is None:
            try:
                response = urllib2.urlopen(pph2_url, params)
                doc = response.read()
            except (socket.timeout, IOError): pass
        soup = BeautifulSoup(doc)
        status_td = soup.find('td', text=re.compile(r'^Batch \d+:'))
        if status_td is None:
            # We might be done, make sure this page is not an error page
            if soup.find('b', text=re.compile(r'^Service Name:')): break
            else:
                tries += 1
                if tries >= max_tries:
                    raise RemoteException('PolyPhen won\'t let us check the status right now.')
                spin(15)
                continue
        pos_td = status_td.parent.parent.findAll('td')[1]
        try: pos = int(pos_td.string)
        except ValueError: pos = 0
        shortened = re.sub(r'^Batch \d+:\s+', '', str(status_td))
        this_step = steps.index(shortened)
        if curr_step != this_step:
            if curr_step is not -1:
                write_status((wait_msg + done_msg) % steps[curr_step], True)
            curr_step += 1
            while curr_step < this_step: # Write out steps that were completed between refreshes.
                write_status((wait_msg + done_msg) % steps[curr_step])
                curr_step += 1
            maxpos = pos
        write_status(wait_msg % shortened, maxpos - pos, maxpos)
        spin(15)
    if curr_step != -1: write_status((wait_msg + done_msg) % steps[curr_step], True)
    else: curr_step = 0
    while curr_step < len(steps): # Write out steps that were completed before last refresh.
        write_status((wait_msg + done_msg) % steps[curr_step])
        curr_step += 1
    result_url = pph2_result_url % sid
    while True:
        error = False
        try:
            write_status("Waiting for PolyPhen2 results => Waiting for download", True)
            response = urllib2.urlopen(result_url)
            result = response.read()
            if result: break
        except (socket.timeout, IOError): spin(15)
    if error: raise RemoteException(error.split("\n")[0])
    write_status(True)
    return result

def update_db_with_results(results_file):
    """ Inserts the contents of the PolyPhen2 result file into our SQLite3 database. """
    rows = results_file.split('\n')
    for row in rows[1:]:
        cells = map(lambda x: x.strip(), row.split('\t'))
        if len(cells) < 17: continue
        tup = (cells[6], cells[7], cells[8], cells[9], cells[16])
        c.execute('INSERT INTO t (accid, seqpos, aa1, aa2, pph2_score) VALUES (?,?,?,?,?)', tup)

def update_db_with_seqlen_and_gene():
    """ For every row in our SQLite3 database, fetch the gene name and sequence length
        from UniProtKB. """
    i = 0
    tot = c.execute('SELECT COUNT(*) FROM t WHERE pph2_score IS NOT NULL').fetchone()[0]
    c.execute('SELECT id, accid FROM t WHERE pph2_score IS NOT NULL')
    for row in c:
        i += 1
        if row is None: continue
        seqlen, gene = get_seqlen_and_gene(str(row[1]))
        write_status("Fetching gene names & protein lengths from UniProtKB", i, tot)
        d.execute('UPDATE t SET seqlen=?, gene=? WHERE id=?', (seqlen, gene, row[0]))
    write_status("Fetching gene names & protein lengths from UniProtKB => Done.\n", True)

def print_results():
    """ Calculates the SUM(pph2_score)/seqlen*1000 per gene, sorts the table,
        and prints it as TSV text to standard out. """
    c.execute("""SELECT gene, (SUM(pph2_score)/seqlen*1000) AS score 
        FROM t WHERE pph2_score IS NOT NULL 
        GROUP BY gene ORDER BY score DESC""")
    csvout = csv.writer(sys.stdout, delimiter='\t')
    for row in c: csvout.writerow(row)

class LocalException(Exception):
    pass

class RemoteException(Exception):
    pass

if __name__=='__main__':  
    import getopt
    help = False
    sid = None
    output = None
    long_args = ["help", "quiet", "sid=", "output="]
    try:
        try:
            opts, args = getopt.getopt(sys.argv[1:], 'ho:qs:', long_args)
        except getopt.GetoptError, err:
            raise LocalException(("Error: %s\n" + USAGE) % (str(err), sys.argv[0]))
        for o, a in opts:
            if o in ("-h", "--help"):
                help = True
            elif o in ("-q", "--quiet"):
                quiet = True
            elif o in ("-o", "--output"):
                output = a
            elif o in ("-s", "--sid"):
                sid = a
            else:
                raise LocalException(("Error: unhandled option: %s\n" + USAGE) % (o, sys.argv[0]))
        if help or (len(args) == 0 and sid is None):
            print USAGE % sys.argv[0]
            sys.exit(1)
        if output is not None:
            sys.stdout = open(output, 'w')
    
        mut_lists = []
        for arg in args:
            if arg=='-': mut_lists.append(sys.stdin.read())
            else: mut_lists.append(open(arg).read())
    
        setup_db()
        if sid is None:
            sid = submit_to_polyphen2("\n".join(mut_lists))
            write_status("Received SID for GGI => %s\n" % sid)
            write_status("Track this job at: %s" % (pph2_track_url % sid))
        results = poll_for_polyphen2_results(sid)
        update_db_with_results(results)
        update_db_with_seqlen_and_gene()
        print_results()
    except (RemoteException, LocalException) as e:
        sys.stderr.write(str(e) + "\n")