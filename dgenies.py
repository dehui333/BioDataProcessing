#!/usr/bin/env python

import argparse
import os
from selenium.common.exceptions import ElementClickInterceptedException
from selenium.common.exceptions import TimeoutException
from selenium.common.exceptions import NoSuchElementException
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
import subprocess
import sys

def startup(port_number, log_path):
    dgenies_log_handle = open(log_path, 'w')
    dgenies_proc = subprocess.Popen(['dgenies', 'run', '-m', 'standalone', '-p', str(port_number), '--no-browser'], stdout=dgenies_log_handle, stderr=dgenies_log_handle)
    return dgenies_log_handle, dgenies_proc

# separate open close browser from plotting?
def plot(port_number, target_path, query_path, output_dir, short_timeout, long_timeout, check_interval):

    target_path = os.path.abspath(target_path)
    query_path = os.path.abspath(query_path)
    output_dir = os.path.abspath(output_dir)
    
    options = Options()
    options.headless = True

    options.set_preference('browser.download.folderList', 2) # custom location
    options.set_preference('browser.download.manager.showWhenStarting', False)
    options.set_preference('browser.download.dir', output_dir)
    options.set_preference('browser.helperApps.neverAsk.saveToDisk', 'text/plain,image/png,text/html')

    driver = webdriver.Firefox(options=options)


    driver.get(f'http://127.0.0.1:{port_number}/')
    
    run_link = driver.find_element(By.LINK_TEXT, 'Run')
    run_link.click()

    # wait object
    short_wait = WebDriverWait(driver, short_timeout)
    check_wait = WebDriverWait(driver, check_interval)

    print('[D-Genies] Moving to run page...', file=sys.stderr)
    

    # wait till page loads
    short_wait.until(EC.presence_of_element_located((By.ID, 'target')))

    
    print('[D-Genies] Filling in forms...', file=sys.stderr)

    # fill in query path
    query_upload = driver.find_element(By.NAME, 'file-query')
    query_upload.send_keys(query_path)

    # fill in target path
    target_upload = driver.find_element(By.NAME, 'file-target')
    target_upload.send_keys(target_path)

    # click submit button
    submit_button = driver.find_element(By.ID, 'submit')
    submit_button.click()

    print('[D-Genies] Waiting for plotting result...', file=sys.stderr)
    # wait for result and go to page

    time_waited = 0
    result_link = None
    while time_waited < long_timeout:
        try:
            result_link = check_wait.until(EC.element_to_be_clickable((By.LINK_TEXT, 'click here')))
            break
        except TimeoutException:
            status = None
            try:
                status = driver.find_element(By.XPATH, "//div[@class='status-body']/p").text
            except:
                print('[D-Genies] Failed to retrieve status!', file=sys.stderr)
            if status:
                print('[D-Genies] Status: ', status, '\n', file=sys.stderr)
            time_waited += check_interval
    
    if result_link == None:
        print('[D-Genies] Timeout while waiting for plotting result!', file=sys.stderr)
        dgenies_log_handle.close()
        driver.close()
        return
    result_link.click()

    #wait
    short_wait.until(EC.presence_of_element_located((By.ID, 'export')))

    print('[D-Genies] Obtaining results...', file=sys.stderr)
    # Select options
    select = Select(driver.find_element(By.XPATH, "//form[@id='export']/select[1]"))
    select.select_by_value('2')
    select.select_by_value('10')
    select.select_by_value('3')
    select.select_by_value('5')
    try: 
        select.select_by_value('6')
    except ElementClickInterceptedException:
        #print('[D-Genies] There is no unmatched queries.', file=sys.stderr)
        pass
    except:
        print('[D-Genies] Something went wrong with retrieving unmatched queries.', file=sys.stderr)

    try: 
        select.select_by_value('7')
    except ElementClickInterceptedException:
        #print('[D-Genies] There is no unmatched targets.', file=sys.stderr)
        pass
    except:
        print('[D-Genies] Something went wrong with retrieving unmatched targets.', file=sys.stderr)

    driver.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot with D-GENIES.')
    parser.add_argument('-i', '--query', type=str, help='Path to query fasta.')
    parser.add_argument('-r', '--target', type=str, help='Path to target fasta.')
    parser.add_argument('-o', '--output', type=str, help='Output directory.')

    args = parser.parse_args()


    plot(5000, args.target, args.query, args.output, 3, 10, 2)
