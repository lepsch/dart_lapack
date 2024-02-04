      import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/matrix.dart';

void dlacpy(final String UPLO, final int M, final int N, final Matrix2d<double> A, final int LDA, final Matrix2d<double> B, final int LDB, ) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // double             A[ LDA, * ], B[ LDB, * ];

      if ( lsame( UPLO, 'U' ) ) {
         for (var J = 1; J <= N; J++) {
            for (var I = 1; I <= min( J, M ); I++) {
               B[I, J] = A[ I, J ];
            }
         }
      } else if ( lsame( UPLO, 'L' ) ) {
         for (var J = 1; J <= N; J++) {
            for (var I = J; I <= M; I++) {
               B[I, J] = A[ I, J ];
            }
         }
      } else {
         for (var J = 1; J <= N; J++) {
            for (var I = 1; I <= M; I++) {
               B[I, J] = A[ I, J ];
            }
         }
      }
      }
