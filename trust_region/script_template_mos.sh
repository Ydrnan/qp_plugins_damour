#!/bin/bash

your_file=

your_C_PROVIDER=
your_H_PROVIDER=
your_g_PROVIDER=
your_cc_PROVIDER=

sed "s/C_PROVIDER/$your_C_PROVIDER/g" trust_region_template_mos.txt > $your_file
sed -i "s/H_PROVIDER/$your_H_PROVIDER/g" $your_file
sed -i "s/g_PROVIDER/$your_g_PROVIDER/g" $your_file
sed -i "s/cc_PROVIDER/$your_cc_PROVIDER/g" $your_file
