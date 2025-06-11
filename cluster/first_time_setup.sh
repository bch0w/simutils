# First time setup on a new cluster, a few standard things that I always do
# written when setting up on Shelob 6/10/25. Some of these instructions are just
# reference and need to be modified to run

UNAME="bhchow"
SNAME="shelob"

# RUN BELOW COMMANDS FROM PERSONAL MACHINE
# ========================================

# Sets up passwordless SSH access
ssh-copy-id ${UNAME}@${SNAME}



# RUN BELOW COMMANDS FROM THE SERVER
# ========================================

# Set up ssh-key
ssh-keygen -t ed25519

# Set up Git account
git config --global user.name bch0w
git config --global user.email "bhchow@alaska.edu"

# Manual Tasks- 
# 1. Copy over personalized bashrc
# 2. Add SSH key to GitHub (https://github.com/settings/ssh/new)
