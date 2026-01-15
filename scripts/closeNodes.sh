#!/bin/bash

WHERETOMOUNT="$HOME/hosts"

# Color codes
GREEN="\033[1;32m"
RED="\033[1;31m"
YELLOW="\033[1;33m"
NC="\033[0m" # No color

echo -e "${YELLOW}Unmounting nodes...${NC}"
sudo umount -l ${WHERETOMOUNT}/nodes
if [ $? -eq 0 ]; then
    echo -e "${GREEN}Nodes data unmounted successfully!${NC}"
else
    echo -e "${RED}Nodes data unmount failed!${NC}"
fi

echo -e "${YELLOW}Unmounting cluster...${NC}"
sudo umount -l ${WHERETOMOUNT}/hpc
if [ $? -eq 0 ]; then
    echo -e "${GREEN}Cluster data unmounted successfully!${NC}"
else
    echo -e "${RED}Cluster data unmount failed!${NC}"
fi

