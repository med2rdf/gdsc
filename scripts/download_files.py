import requests
import os
from urllib.parse import urlparse
from ftplib import FTP


def download(url, file_name, dir_name):
    res = requests.get(url)
    file_path = os.path.join(dir_name, file_name)
    if os.path.exists(file_path):
        return
    print('download from: {}'.format(url))
    with open(file_path, 'w') as f:
        f.write(res.text)
    return


def download_via_ftp(url, file_name, dir_name):
    print('download from: {}'.format(url))
    file_path = os.path.join(dir_name, file_name)
    if os.path.exists(file_path):
        return
    host = urlparse(url).netloc
    ftp_file_path = urlparse(url).path
    ftp = FTP(host)
    ftp.login()
    with open(file_path, "w") as f:
        ftp.retrlines('RETR ' + ftp_file_path, f.write)
    ftp.quit()
    return


save_dir = 'data/raw'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

download('https://www.cancerrxgene.org/downloads/download/ic', 'PANCANCER_IC.csv', save_dir)
download('https://www.cancerrxgene.org/downloads/download/anova', 'PANCANCER_ANOVA.csv', save_dir)
download_via_ftp('ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt', 'cellosaurus.txt', save_dir)
