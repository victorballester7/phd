#!/bin/bash

NODE="typhoon" # the node actually does not play a role because the databases are all linked, so it's the same for all nodes: typhoon, hurricane, spitfire, blackfriars.
CLUSTER="hpc"
USER="vb824"

FOLDERTOMOUNT="Desktop/"
WHERETOMOUNT="$HOME/hosts"

# Color codes
GREEN="\033[1;32m"
RED="\033[1;31m"
YELLOW="\033[1;33m"
NC="\033[0m" # No color

echo -e "${YELLOW}Mounting nodes...${NC}"
mkdir -p ${WHERETOMOUNT}/nodes
sshfs -o allow_other -o kernel_cache -o auto_cache -o reconnect -o compression=no -o cache_timeout=600 -o ServerAliveInterval=15 ${USER}@${NODE}:/home/${USER}/${FOLDERTOMOUNT} ${WHERETOMOUNT}/nodes
if [ $? -eq 0 ]; then
    echo -e "${GREEN}Nodes data mounted successfully!${NC}"
else
    echo -e "${RED}Nodes data mount failed!${NC}"
fi

echo -e "${YELLOW}Mounting cluster...${NC}"
mkdir -p ${WHERETOMOUNT}/hpc
sshfs -o allow_other -o kernel_cache -o auto_cache -o reconnect -o compression=no -o cache_timeout=600 -o ServerAliveInterval=15 ${USER}@${CLUSTER}:/rds/general/user/${USER}/home/${FOLDERTOMOUNT} ${WHERETOMOUNT}/hpc
if [ $? -eq 0 ]; then
    echo -e "${GREEN}Cluster data mounted successfully!${NC}"
else
    echo -e "${RED}Cluster data mount failed!${NC}"
fi

