#ifndef  IBLAS_H
#define  IBLAS_H
#include "conv/usconv.h"
#include "uscr/uscr.h"
#include "utils.h"
#ifdef SUPPORT_MIC
	#include "mic/mic_mbv.h"
#else
	#ifdef SUPPORT_GPU
		#include "gpu/gpu_mbv.h"
	#else
		#include "host/host_mbv.h"
	#endif
#endif
#endif