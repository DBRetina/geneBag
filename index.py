import kProcessor as kp
import sys

kSize = 25

if len(sys.argv) < 2:
    sys.exit("run: python index.py <fasta_file> <kSize:default(25)>")

elif len(sys.argv) == 2:
    fasta_file = sys.argv[1]

elif len(sys.argv) == 3:
    fasta_file = sys.argv[1]
    kSize = int(sys.argv[2])

chunkSize = 5000

KF = kp.kDataFramePHMAP(kp.KMERS, kp.nonCanonicalInteger_Hasher, {"kSize":kSize})
cKF = kp.index(KF, fasta_file, chunkSize, fasta_file + ".names")
cKF.save("idx_" + fasta_file.replace(".fa",""))
