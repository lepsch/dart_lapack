import 'dart:math';

import 'package:lapack/src/blas/lsame.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/ilaenv.dart';
import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/xerbla.dart';

      void dsyswapr(UPLO, N, A, LDA, I1, I2) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String           UPLO;
      int              I1, I2, LDA, N;
      double           A( LDA, * );

// =====================================================================

      bool               UPPER;
      double             TMP;

      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DSWAP

      UPPER = lsame( UPLO, 'U' );
      if (UPPER) {

          // UPPER
          // first swap
          //  - swap column I1 and I2 from I1 to I1-1
         dswap(I1-1, A(1,I1), 1, A(1,I2), 1 );

           // second swap :
           // - swap A(I1,I1) and A(I2,I2)
           // - swap row I1 from I1+1 to I2-1 with col I2 from I1+1 to I2-1
         TMP=A(I1,I1);
         A(I1,I1)=A(I2,I2);
         A(I2,I2)=TMP;

         dswap(I2-I1-1, A(I1,I1+1), LDA, A(I1+1,I2), 1 );

           // third swap
           // - swap row I1 and I2 from I2+1 to N
         if (I2 < N) dswap( N-I2, A(I1,I2+1), LDA, A(I2,I2+1), LDA );

        } else {

          // LOWER
          // first swap
          //  - swap row I1 and I2 from I1 to I1-1
         dswap(I1-1, A(I1,1), LDA, A(I2,1), LDA );

          // second swap :
          //  - swap A(I1,I1) and A(I2,I2)
          //  - swap col I1 from I1+1 to I2-1 with row I2 from I1+1 to I2-1
          TMP=A(I1,I1);
          A(I1,I1)=A(I2,I2);
          A(I2,I2)=TMP;

          dswap(I2-I1-1, A(I1+1,I1), 1, A(I2,I1+1), LDA );

          // third swap
          //  - swap col I1 and I2 from I2+1 to N
         if (I2 < N) dswap( N-I2, A(I2+1,I1), 1, A(I2+1,I2), 1 );

      }
      }
