#!/usr/bin/env python

import sys, os, bvcp_core

read_data = (sys.argv[3], sys.argv[4], sys.argv[2], sys.argv[5], sys.argv[6], os.path.dirname(sys.argv[7]))
bvcp_core.run_smalt(read_data, sys.argv[1], os.path.join(sys.argv[1],"intermediate_files"), sys.argv[10], sys.argv[8], sys.argv[9], sys.argv[7])
