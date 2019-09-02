import requests
import os
from urllib.parse import urlparse
from ftplib import FTP


def download(url, save_dir):
    print('download from: {}'.format(url))
    res = requests.get(url)
    dispositoin = res.headers['content-disposition']
    file_name = dispositoin.split('filename=')[-1].strip('"')
    with open(os.path.join(save_dir, file_name), 'w') as f:
        f.write(res.text)
    return


def download_via_ftp(url, save_dir):
    print('download from: {}'.format(url))
    host = urlparse(url).netloc
    ftp_file_path = urlparse(url).path
    file_name = os.path.split(ftp_file_path)[-1]

    ftp = FTP(host)
    ftp.login()
    with open(os.path.join(save_dir, file_name), "w") as f:
        ftp.retrlines('RETR ' + ftp_file_path, f.write)
    ftp.quit()
    return


save_dir = 'data/raw'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

download('https://www.cancerrxgene.org/downloads/download/ic', save_dir)
download('https://www.cancerrxgene.org/downloads/download/anova', save_dir)
download_via_ftp('ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt', save_dir)

