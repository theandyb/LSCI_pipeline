"""
Store reference genome as sqlite db
"""

import sqlite3
from Bio import SeqIO

# Function for the creation of the database
def create_database(fasta_file, db_file):
  conn = sqlite3.connect(db_file)
  c = conn.cursor()
  
  c.execute('''
  CREATE TABLE IF NOT EXISTS sequences (
    id TEXT PRIMARY KEY,
    sequence TEXT
  )
  ''')
  
  c.execute('''CREATE INDEX IF NOT EXISTS idx_name ON sequences (id)''')
  
  for record in SeqIO.parse(fasta_file, "fasta"):
    c.execute('INSERT INTO sequences (id, sequence) VALUES (?, ?)', (record.id, str(record.seq)))
  
  conn.commit()
  conn.close()

fasta_file="data/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
db_file="data/reference/ref_db.db"

create_database(fasta_file, db_file)
