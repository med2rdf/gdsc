import requests
import os
from urllib.parse import urlparse
from ftplib import FTP
from progressbar import ProgressBar, UnknownLength


def download(url, dir_name, file_name):
    download_path = os.path.join(dir_name, file_name)
    if os.path.exists(download_path):
        return
    scheme = urlparse(url).scheme
    if scheme in ['http', 'https']:
        download_with_progress_bar(url, download_path)
    elif scheme == 'ftp':
        download_via_ftp_with_progress_bar(url, download_path)
    return


def download_with_progress_bar(url, download_path):
    print('download from: {}'.format(url))
    file_size = int(requests.head(url).headers.get('content-length', 0))
    if file_size == 0:
        file_size = UnknownLength
    res = requests.get(url, stream=True)

    pbar = ProgressBar(maxval=file_size)
    pbar.start()
    with open(download_path, 'wb') as f:
        for chunk in res.iter_content(chunk_size=1024):
            f.write(chunk)
            pbar.update(len(chunk))
    pbar.finish()
    return


def download_via_ftp_with_progress_bar(url, download_path):
    print('download from: {}'.format(url))
    host = urlparse(url).netloc
    ftp_file_path = urlparse(url).path
    ftp = FTP(host)
    ftp.login()
    try:
        file_size = ftp.size(ftp_file_path)
    except:
        file_size = UnknownLength
    pbar = ProgressBar(maxval=file_size)
    pbar.start()
    with open(download_path, 'w') as f:
        ftp.retrlines('RETR ' + ftp_file_path, lambda data: update_pbar(data+'\n', f, pbar))
    pbar.finish()
    ftp.quit()
    return


def update_pbar(data, f, pbar):
    f.write(data)
    pbar += len(data)
    return


save_dir = 'data/raw'
if not os.path.isdir(save_dir):
    os.makedirs(save_dir)

download('https://www.cancerrxgene.org/downloads/download/ic', save_dir, 'PANCANCER_IC.csv')
download('https://www.cancerrxgene.org/downloads/download/anova', save_dir, 'PANCANCER_ANOVA.csv')
download('ftp://ftp.expasy.org/databases/cellosaurus/cellosaurus.txt', save_dir, 'cellosaurus.txt')
