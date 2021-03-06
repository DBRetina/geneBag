import kProcessor as kp
import sys

kSize = 25
fasta_file = str()
chunkSize = 5000

if len(sys.argv) < 2:
    sys.exit("run: python index.py <fasta_file> <kSize:default(25)>")

elif len(sys.argv) == 2:
    fasta_file = sys.argv[1]

elif len(sys.argv) == 3:
    fasta_file = sys.argv[1]
    kSize = int(sys.argv[2])

names_file = fasta_file + ".names"

# want it fast but more memory cost?
KF = kp.kDataFramePHMAP(kp.KMERS, kp.nonCanonicalInteger_Hasher, {"kSize":kSize})

# a bit slower but memory effiecent?
#KF = kp.kDataFrameMQF(kp.KMERS, kp.nonCanonicalInteger_Hasher, {"kSize":kSize})

# MQF with high intial memory 
#KF = kp.kDataFrameMQF(kSize, 34, kp.nonCanonicalInteger_Hasher)

cKF = kp.index(KF, fasta_file, chunkSize, names_file)
cKF.save("idx_" + fasta_file.replace(".fa",""))
