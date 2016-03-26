vpath %.h include
vpath %.h include/host
vpath %.h include/uscr
vpath %.h include/conv

vpath %.c src
vpath %.c src/conv
vpath %.c src/host
vpath %.c src/sbv
vpath %.c src/uscr

CFLAGS = -g -fopenmp -lm

objects=sp_io.o host_usmv.o host_rmbv_bco.o host_rmbv_bdi.o host_rmbv_bsc.o host_rmbv_bsr.o host_rmbv_coo.o \
host_rmbv_csc.o host_rmbv_csr.o host_rmbv_vbr.o host_rmbv_dia.o link.o dense.o properties.o \
comm_tools.o inserting.o host_lmbv_bco.o host_lmbv_bdi.o host_lmbv_bsc.o host_lmbv_bsr.o \
host_lmbv_coo.o host_lmbv_csc.o host_lmbv_csr.o host_lmbv_vbr.o host_lmbv_dia.o host_usmm.o usds.o \
uscr_coo.o uscr_csc.o uscr_csr.o uscr_bsc.o uscr_bsr.o uscr_bco.o \
uscr_dia.o uscr_bdi.o uscr_vbr.o ins_routiner.o uscr_begin.o  \
uscr_block_begin.o uscr_variable_block_begin.o uscr_end.o uscr_insert_entry.o \
uscr_insert_block.o uscr_insert_row.o uscr_insert_col.o uscr_insert_clique.o \
uscr_insert_entries.o conv_tools.o usconv_bco2bdi.o usconv_bco2bsc.o usconv_bco2bsr.o usconv_bco2coo.o \
usconv_bdi2bco.o usconv_bsc2bco.o usconv_bsr2bco.o usconv_coo2csr.o usconv_coo2csc.o \
usconv_coo2dia.o usconv_coo2bco.o usconv_csc2coo.o usconv_csr2bco.o usconv_csr2coo.o usconv_dia2coo.o usaxpy.o usdot.o \
usga.o usgz.o ussc.o mmio.o

lib:$(objects)
	ar crv libspblas.a $(objects)

sp_io.o:sp_io.h
host_usmv.o:host_USMV.h
host_rmbv_bco.o:host_rmbv_bco.h dense.h 
host_rmbv_bdi.o:host_rmbv_bdi.h dense.h 
host_rmbv_bsc.o:host_rmbv_bsc.h dense.h 
host_rmbv_bsr.o:host_rmbv_bsr.h dense.h 
host_rmbv_coo.o:host_rmbv_coo.h link.h 
host_rmbv_csc.o:host_rmbv_csc.h 
host_rmbv_csr.o:host_rmbv_csr.h 
host_rmbv_vbr.o:host_rmbv_vbr.h dense.h 
host_rmbv_dia.o:host_rmbv_dia.h
host_lmbv_bco.o:host_lmbv_bco.h dense.h 
host_lmbv_bdi.o:host_lmbv_bdi.h dense.h 
host_lmbv_bsc.o:host_lmbv_bsc.h dense.h 
host_lmbv_bsr.o:host_lmbv_bsr.h dense.h 
host_lmbv_coo.o:host_lmbv_coo.h
host_lmbv_csc.o:host_lmbv_csc.h 
host_lmbv_csr.o:host_lmbv_csr.h
host_lmbv_vbr.o:host_lmbv_vbr.h dense.h 
host_lmbv_dia.o:host_lmbv_dia.h

host_usmm.o:host_USMM.h
usds.o:usds.h

uscr_coo.o:uscr_coo.h
uscr_csc.o:uscr_csc.h
uscr_csr.o:uscr_csr.h
uscr_bsc.o:uscr_bsc.h
uscr_bsr.o:uscr_bsr.h
uscr_bco.o:uscr_bco.h
uscr_dia.o:uscr_dia.h types.h comm_tools.h properties.h usds.h
uscr_bdi.o:uscr_bdi.h
uscr_vbr.o:uscr_vbr.h
ins_routiner.o:INS_ROUTINER.h
uscr_begin.o:USCR_BEGIN.h
uscr_block_begin.o:USCR_BEGIN.h
uscr_variable_block_begin.o:USCR_BEGIN.h
uscr_end.o:USCR_END.h
uscr_insert_entry.o:USCR_INSERT_ENTRY.h
uscr_insert_block.o:USCR_INSERT_BLOCK.h USCR_INSERT_ENTRY.h
uscr_insert_row.o:USCR_INSERT_ENTRY.h USCR_INSERT_ROW.h
uscr_insert_col.o:USCR_INSERT_ENTRY.h USCR_INSERT_COL.h
uscr_insert_clique.o:USCR_INSERT_ENTRY.h USCR_INSERT_CLIQUE.h
uscr_insert_entries.o:USCR_INSERT_ENTRY.h USCR_INSERT_ENTRIES.h
conv_tools.o:conv_tools.h blas_sparse_namedconstants.h

usconv_bco2bdi.o:properties.h conv_tools.h usconv.h link.h
usconv_bco2bsc.o:properties.h conv_tools.h usconv.h link.h
usconv_bco2bsr.o:properties.h conv_tools.h usconv.h link.h
usconv_bco2coo.o:properties.h conv_tools.h usconv.h link.h
usconv_bdi2bco.o:properties.h conv_tools.h usconv.h link.h
usconv_bsc2bco.o:properties.h conv_tools.h usconv.h link.h
usconv_bsr2bco.o:properties.h conv_tools.h usconv.h link.h
usconv_coo2csr.o:properties.h conv_tools.h usconv.h link.h
usconv_coo2csc.o:properties.h conv_tools.h usconv.h link.h
usconv_coo2dia.o:properties.h conv_tools.h usconv.h link.h
usconv_coo2bco.o:properties.h conv_tools.h usconv.h link.h
usconv_csc2coo.o:properties.h conv_tools.h usconv.h link.h
usconv_csr2bco.o:properties.h conv_tools.h usconv.h link.h
usconv_csr2coo.o:properties.h conv_tools.h usconv.h link.h
usconv_dia2coo.o:properties.h conv_tools.h usconv.h link.h

link.o:link.h
dense.o:dense.h blas_sparse_namedconstants.h
properties.o:properties.h
comm_tools.o:comm_tools.h INSERTING.h properties.h types.h link.h
inserting.o:INSERTING.h comm_tools.h
usaxpy.o:USAXPY.h
usdot.o:USDOT.h comm_tools.h
usga.o:USGA.h
usgz.o:USGA.h USGZ.h
ussc.o:USSC.h
mmio.o:mmio.h
#hash.o:hash.h

# test.o:USMV.h USMM.h USCR_BEGIN.h USCR_END.h USCR_INSERT_ENTRY.h USCR_INSERT_BLOCK.h usds.h comm_tools.h uscr_coo.h \
# 	uscr_csr.h uscr_bsc.h uscr_bsr.h uscr_bco.h uscr_dia.h uscr_bdi.h uscr_vbr.h mbv.h blas_enum.h usconv.h link.h
# 	gcc -c -o test.o mytest.c
# gcc -shared -fPCI -o libspblas.so $(objects)
# lib:$(objects)
# 	ar crv libspblas.a $(objects)
mytest:
	gcc -g -fopenmp -lm -o mytest mytest.c libspblas.a
test:
	gcc -g -fopenmp -lm -o test test.c libspblas.a
temp:
	gcc -g -fopenmp -lm -o temp_property temp_property.c libspblas.a
feature:
	gcc -g -fopenmp -lm -o feature feature.c libspblas.a
	

.PHONE:clean.o clean clean.a clean.so
clean:
	rm mytest test feature temp_property $(objects) libspblas.a
clean.o:
	rm $(objects)