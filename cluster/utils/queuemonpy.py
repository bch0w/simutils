"""
A daemon to monitor jobs on a SLURM HPC system with notifications via email.
Allows for fire-and-forget job running without having to log onto the cluster
to check job status. Designed for the NeSI system Maui.
"""
import sys
import time
import smtplib
import argparse
import subprocess

def parse_args():
    """
    General argument parser to provide a level of customization
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("jobid", type=str, nargs="?",
                        help="The numerical job identifier, collected as list")
    parser.add_argument("-p", "--password", required=True, type=str, nargs="?")
    parser.add_argument("-a", "--account", type=str, nargs=1,
                        default="chowbr", help="Account name on cluster")
    parser.add_argument("-c", "--cluster", type=str, nargs=1,
                        default="maui,maui_ancil", help="Cluster to monitor")
    parser.add_argument("-u", "--user", type=str, nargs=1,
                        default="queuemonpy@gmail.com",
                        help="Gmail user account for email sending")
    parser.add_argument("-w", "--wait", type=int, nargs=1, default=60,
                        help="Time to wait between sacct queries to cluster in "
                             "minutes")
    parser.add_argument("-i", "--initial", nargs="?", type=int, default=0,
                        help="Send an initial email when tracking starts, "
                             "0 for no, 1 for yes")

    return parser.parse_args()


class Queuemonpy:
    """
    A class that monitors jobs and sends emails
    """
    def __init__(self, check_login=True):
        """
        Simply parse arguments from command line. Check that email login works
        before finishing initialization.
        """
        self.args = parse_args()
        if check_login:
            try:
                if bool(self.args.initial):
                    self.send(subject=f"QMONPY TRACKING: {self.args.jobid}",
                              body="Courtesy email that job tracking has begun")
                else:
                    self.login()
            except Exception as e:
                print(f"Problem signing into email account "
                      f"{self.args.user}:\n{e}")
                sys.exit(-1)

    def login(self):
        """
        Login to the email client, used for testing the email client before
        beginning job tracking
        """
        server = smtplib.SMTP_SSL("smtp.gmail.com", 465)
        server.ehlo()
        server.login(self.args.user, self.args.password)
        server.close()

    def send(self, subject, body):
        """
        Send an email to and from the given user with the given subject and body
        Requires that less than secure apps are allowed access to the gmail 
        account. Recommended that a non-critical email is used.

        :type subject: str
        :param subject: subject or header line of the email
        :type body: str
        :param body: body of the email
        """
        email_text = (f"From {self.args.user}\n"
                      f"To: {self.args.user}\n"
                      f"Subject: {subject}\n"
                      f"\n"
                      f"{body}")

        server = smtplib.SMTP_SSL("smtp.gmail.com", 465)
        server.ehlo()
        server.login(self.args.user, self.args.password)
        server.sendmail(self.args.user, self.args.user, email_text)
        server.close()

    def check(self, jobid):
        """
        Check status on a single job using subprocess and the 'sacct' SLURM

        :type jobid: str
        :param jobid: job id to check using sacct
        :rtype: str or None
        :return: str corresponding to the status, if no matching job found,
            returns None for status
        """
        check = (f"sacct -nL --cluster={self.args.cluster} -o jobid,state " \
                 f"-j {jobid}")
        stdout = subprocess.run(check, capture_output=True, shell=True,
                                text=True).stdout
        final_status = None
        for job in stdout.strip().split("\n"):
            jobid_temp, status = job.split()
            if jobid_temp == jobid:
                # Sometimes there are two matching jobids, take the latest one
                final_status = status
            else:
                # Any jobs that do not match the job id exactly are auxiliary
                continue
        return final_status

    def sacct(self, jobid):
        """
        Return the full sacct output from a job id, to be put into the body of
        the email sent by the daemon

        :type jobid: str
        :param jobid: job id to check using sacct
        :rtype: str or None
        :return: str corresponding to the status, if no matching job found,
            returns None for status
        """
        check = f"sacct --clusters={self.args.cluster} -j {jobid}"
        return subprocess.run(check, capture_output=True, shell=True,
                              text=True).stdout

    def monitor(self):
        """
        Monitor ongoing jobs by repeatedly checking the status of jobsv until
        status other than 'RUNNING' is returned.
        """
        jobids = [_ for _ in self.args.jobid.split(",")]
        while jobids:
            # [:] ensures that we can remove from an iterating list
            for jobid in jobids[:]:
                status = self.check(jobid)
                if status == "RUNNING":
                    continue
                else:
                    self.send(subject=f"{jobid} {status}",
                              body=self.sacct(jobid))
                    jobids.remove(jobid)

            time.sleep(self.args.wait)

if __name__ == "__main__":
    qmon = Queuemonpy()
    qmon.monitor()
