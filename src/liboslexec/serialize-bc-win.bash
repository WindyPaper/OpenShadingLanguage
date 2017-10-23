echo "#include <cstddef>" > $2
echo "unsigned char osl_llvm_compiled_ops_block[] = {" >> $2
E:/cygwin64/bin/hexdump -v -e '"" /1 "0x%02x" ",\n"' $1 >> $2
echo "0x00 };" >> $2
echo "size_t osl_llvm_compiled_ops_size = sizeof(osl_llvm_compiled_ops_block)-1;" >> $2
read -t5 -n1 -r -p 'Press any key in the next five seconds...' key