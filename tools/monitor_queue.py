"""
Convenience function to check if jobs are still running so that I don't have 
to run the squeue function over and over again. Based off Seisflows job_status
"""
import sys
import time
import smtplib
import subprocess


def monitor(wait_time_s=60, cluster_names="maui,maui_ancil"):
    """
    Monitor current jobs on cluster using slurm command squeue

    :type wait_time_s: int
    :param wait_time_s: time to wait before rechecking status
    :type cluster_names: str
    :param cluster_names: list of cluster names separated by commas
    """
    check_command = 'squeue --user=chowbr --clusters={}'
    print_job = "{c} has {j} jobs; {n}; time_left={t}"

    # Keep loop open until all jobs have finished
    while True:
        clusters = {}

        for cluster_name in cluster_names.split(','):
            output = subprocess.run(check_command.format(cluster_name).split(), 
                                    capture_output=True, text=True)
            lines = output.stdout.strip().split('\n')

            jobids, names, time_lefts = [],[],[]
            for line in lines:
                if 'CLUSTER' in line:
                    cluster = line.split(':')[1].strip()
                elif 'JOBID' in line:
                    continue
                else:
                    jobid, _, _, name, _, _, _, _, time_left, _, _ = \
                                                                    line.split()
                    time_lefts.append(time_left)
                    jobids.append(jobid)
                    names.append(name)

            if jobids:
                clusters[cluster] = {"njobids":len(jobids),
                                     "jobids": jobids,
                                     "names":list(set(names)),
                                     "maxtime":max(time_lefts)
                                     }
       
        # parse the dictionaries created 
        still_jobs = False
        for ckey in clusters.keys():
            if 'jobids' in clusters[ckey].keys():
                still_jobs = True
                # if job number stays the same, don't print
                try:
                    if clusters[ckey]['jobids'] == old_clusters[ckey]['jobids']:
                        continue
                except (NameError, KeyError):
                    pass

                # print some relevant information for the user to see job status
                print(print_job.format(c=ckey, j=clusters[ckey]['njobids'], 
                                       n=clusters[ckey]['names'], 
                                       t=clusters[ckey]['maxtime']
                                       ) 
                      )

        
        # check to see if loop should be run again
        if still_jobs:
            old_clusters = clusters
            time.sleep(wait_time_s)   
        else:
            print("No Jobs on {}".format(cluster_names.split(',')))
            return True


def mailer(subject, body):
    """
    Send an email using Gmail

    :type subject: str
    :param subject: subject of the email
    :type body: str
    :param body: body of the email
    """
    gmail_user = "bch0w.daemon@gmail.com"
    gmail_password = "gkfgdfjukckekcgv"
    
    email_text = f"""\
From: {gmail_user}
To: {gmail_user}
Subject: {subject}
    
{body}
    """
    
    try:
        server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
        server.ehlo()
        server.login(gmail_user, gmail_password)
        server.sendmail(gmail_user, gmail_user, email_text)
        server.close()   
    except:
        print("mail fail")


def monitor_master_job(jobids, mail=True):
    """
    Seisflows will run a master job in the background. When this job ends
    either naturally or through an error, send an email using the mailer 

    :type jobids: str or list of str
    :param jobids: jobids to be checked
    :type mail: bool
    :param mail: mail if job stops running
    """
    jobids = jobids.split(',')
    while True:
        # make sure jobids is a list object
        for jobid in jobids:
            check_command = f'sacct -nL -o jobname,state -j {jobid}'
            output = subprocess.run(check_command.split(), capture_output=True, 
                                    text=True)

            # slurm has auxiliary jobs that run under the same job number,
            # so I named my master jobs with '_master' for easy id
            lines = output.stdout.split()
            for i, line in enumerate(lines):
                # !!!  this needs to be changed, cannot fit _master into jobname
                if "_master" in line:  
                    jobname = line
                    jobstatus = lines[i+1]
                    break
            # job is still running, move on
            if jobstatus == "RUNNING":
                continue
            else:
                if mail:
                    mailer(subject=f"{jobid} {jobstatus}", body="")
                jobids.remove(jobid)
                if not jobids:
                    print("Jobs not running")
                    return
        
        # Wait to check again
        time.sleep(300)

if __name__ == "__main__":
    try:
        jobids = sys.argv[1]
        monitor_master_job(jobids)
    except IndexError:
        job_type = input("[l]ive monitoring or [a]fk?: ")
        if job_type == "l":
            monitor()
        elif job_type == "a":
            jobid = input("jobid(s)?: ")
            monitor_master_job(jobid)
        else:
            print('')
        
