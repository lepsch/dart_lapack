// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

library;

export 'base.dart';
export 'intrinsics.dart';

// Double precision real BLAS 1 routines
export 'src/blas/idamax.dart';
export 'src/blas/dasum.dart';
export 'src/blas/daxpy.dart';
export 'src/blas/dcopy.dart';
export 'src/blas/ddot.dart';
export 'src/blas/dnrm2.dart';
export 'src/blas/drot.dart';
export 'src/blas/drotg.dart';
export 'src/blas/dscal.dart';
export 'src/blas/dsdot.dart';
export 'src/blas/dswap.dart';
export 'src/blas/drotmg.dart';
export 'src/blas/drotm.dart';

// ZBLAS1 - Double precision complex BLAS routines
export 'src/blas/dcabs1.dart';
export 'src/blas/dzasum.dart';
export 'src/blas/dznrm2.dart';
export 'src/blas/izamax.dart';
export 'src/blas/zaxpy.dart';
export 'src/blas/zcopy.dart';
export 'src/blas/zdotc.dart';
export 'src/blas/zdotu.dart';
export 'src/blas/zdscal.dart';
export 'src/blas/zrotg.dart';
export 'src/blas/zscal.dart';
export 'src/blas/zswap.dart';
export 'src/blas/zdrot.dart';

// ALLBLAS - Auxiliary routines for Level 2 and 3 BLAS
export 'src/blas/lsame.dart';
export 'src/blas/xerbla.dart';
export 'src/blas/xerbla_array.dart';

// DBLAS2 - Double precision real BLAS2 routines
export 'src/blas/dgemv.dart';
export 'src/blas/dgbmv.dart';
export 'src/blas/dsymv.dart';
export 'src/blas/dsbmv.dart';
export 'src/blas/dspmv.dart';
export 'src/blas/dtrmv.dart';
export 'src/blas/dtbmv.dart';
export 'src/blas/dtpmv.dart';
export 'src/blas/dtrsv.dart';
export 'src/blas/dtbsv.dart';
export 'src/blas/dtpsv.dart';
export 'src/blas/dger.dart';
export 'src/blas/dsyr.dart';
export 'src/blas/dspr.dart';
export 'src/blas/dsyr2.dart';
export 'src/blas/dspr2.dart';

// ZBLAS2 - Double precision complex BLAS2 routines
export 'src/blas/zgemv.dart';
export 'src/blas/zgbmv.dart';
export 'src/blas/zhemv.dart';
export 'src/blas/zhbmv.dart';
export 'src/blas/zhpmv.dart';
export 'src/blas/ztrmv.dart';
export 'src/blas/ztbmv.dart';
export 'src/blas/ztpmv.dart';
export 'src/blas/ztrsv.dart';
export 'src/blas/ztbsv.dart';
export 'src/blas/ztpsv.dart';
export 'src/blas/zgerc.dart';
export 'src/blas/zgeru.dart';
export 'src/blas/zher.dart';
export 'src/blas/zhpr.dart';
export 'src/blas/zher2.dart';
export 'src/blas/zhpr2.dart';

// DBLAS3 - Double precision real BLAS3 routines
export 'src/blas/dgemm.dart';
export 'src/blas/dsymm.dart';
export 'src/blas/dsyrk.dart';
export 'src/blas/dsyr2k.dart';
export 'src/blas/dtrmm.dart';
export 'src/blas/dtrsm.dart';

// ZBLAS3 - Double precision complex BLAS3 routines
export 'src/blas/zgemm.dart';
export 'src/blas/zsymm.dart';
export 'src/blas/zsyrk.dart';
export 'src/blas/zsyr2k.dart';
export 'src/blas/ztrmm.dart';
export 'src/blas/ztrsm.dart';
export 'src/blas/zhemm.dart';
export 'src/blas/zherk.dart';
export 'src/blas/zher2k.dart';
