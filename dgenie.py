#!/usr/bin/env python

import argparse
import os
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.firefox.options import Options
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import Select
from selenium.webdriver.support.ui import WebDriverWait
import sys


def run():
    pass

def plot(port_number, target_path, query_path, output_dir, timeout_seconds):
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
    wait = WebDriverWait(driver, timeout_seconds)

    print('[D-DEGENIES] Moving to run page...', file=sys.stderr)
    

    # wait till page loads
    wait.until(EC.presence_of_element_located((By.ID, 'target')))

    
    print('[D-DEGENIES] Filling in forms...', file=sys.stderr)

    # fill in query path
    query_upload = driver.find_element(By.NAME, 'file-query')
    query_upload.send_keys(query_path)

    # fill in target path
    target_upload = driver.find_element(By.NAME, 'file-target')
    target_upload.send_keys(target_path)

    # click submit button
    submit_button = driver.find_element(By.ID, 'submit')
    submit_button.click()

    print('[D-DEGENIES] Waiting for plotting result...', file=sys.stderr)
    # wait for result and go to page
    result_link = wait.until(EC.element_to_be_clickable((By.LINK_TEXT, 'click here')))
    result_link.click()

    #wait
    wait.until(EC.presence_of_element_located((By.ID, 'export')))

    print('[D-DEGENIES] Obtaining results...', file=sys.stderr)
    # Select options
    select = Select(driver.find_element(By.XPATH, "//form[@id='export']/select[1]"))
    select.select_by_value('2')
    select.select_by_value('10')
    select.select_by_value('3')
    select.select_by_value('5')
    select.select_by_value('6')
    select.select_by_value('7')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Plot with D-GENIES.')
    parser.add_argument('-i', '--query', type=str, help='Path to query fasta.')
    parser.add_argument('-r', '--target', type=str, help='Path to target fasta.')
    parser.add_argument('-o', '--output', type=str, help='Output directory.')

    args = parser.parse_args()
    if os.path.isdir(args.output):
        print(f'{args.output} already exists!', file=sys.stderr)
        exit(1)
    os.mkdir(args.output)

    plot(5000, os.path.abspath(args.target), os.path.abspath(args.query), os.path.abspath(args.output), 600)
