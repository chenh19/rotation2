#!/bin/bash

# check for existing keys
#ls -al ~/.ssh

# ask for your_email
read -p 'Please enter your GitHub email address: ' email


# create a key if does not exist
ssh-keygen -t ed25519 -C $email
echo 'when asked "Enter a file in which to save the key," press Enter (default file location)'
echo 'input a passphrase (anything you can remember)'

# add SSH key to ssh-agent-
eval `ssh-agent -s`
ssh-add ~/.ssh/ed25519

# get the key
cat ~/.ssh/ed25519.pub
echo 'copy the above key, then navigate to Github account and add the ssh key'
