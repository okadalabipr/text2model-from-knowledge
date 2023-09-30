from ftplib import FTP
import os

if not os.path.exists("data"):
    os.makedirs("data")

# download data from PubTator Central ftp server
url = "ftp.ncbi.nlm.nih.gov"
path = "/pub/lu/PubTatorCentral/PubTatorCentral_BioCXML/"

ftp = FTP(url)
ftp.login()
ftp.cwd(path)

filenames = ftp.nlst()

for filename in filenames:
    with open(os.path.join("data", filename), "wb") as f:
        ftp.retrbinary("RETR " + filename, f.write)
ftp.quit()
