      void zdrvsg2stg(NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, LDB, D, D2, Z, LDZ, AB, BB, AP, BP, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, RESULT, INFO ) {

      // IMPLICIT NONE

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LDZ, LIWORK, LRWORK, NOUNIT, NSIZES, NTYPES, NWORK;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double             D( * ), D2( * ), RESULT( * ), RWORK( * );
      Complex         A( LDA, * ), AB( LDA, * ), AP( * ), B( LDB, * ), BB( LDB, * ), BP( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ZERO, ONE, TEN;
      const              ZERO = 0.0, ONE = 1.0, TEN = 10.0 ;
      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             UPLO;
      int                I, IBTYPE, IBUPLO, IINFO, IJ, IL, IMODE, ITEMP, ITYPE, IU, J, JCOL, JSIZE, JTYPE, KA, KA9, KB, KB9, M, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double             ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU, TEMP1, TEMP2;
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLARND;
      // EXTERNAL lsame, DLAMCH, DLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLAFTS, DLASUM, XERBLA, ZHBGV, ZHBGVD, ZHBGVX, ZHEGV, ZHEGVD, ZHEGVX, ZHPGV, ZHPGVD, ZHPGVX, ZLACPY, ZLASET, ZLATMR, ZLATMS, ZSGT01, ZHEGV_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 9,];
      const KMAGN = [ 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 1, 1, 1,];
      const KMODE = [ 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4,];
      // ..
      // .. Executable Statements ..

      // 1)      Check for errors

      NTESTT = 0;
      INFO = 0;

      BADNN = false;
      NMAX = 0;
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = max( NMAX, NN( J ) );
         if( NN( J ) < 0 ) BADNN = true;
      } // 10

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1;
      } else if ( BADNN ) {
         INFO = -2;
      } else if ( NTYPES < 0 ) {
         INFO = -3;
      } else if ( LDA <= 1 || LDA < NMAX ) {
         INFO = -9;
      } else if ( LDZ <= 1 || LDZ < NMAX ) {
         INFO = -16;
      } else if ( 2*max( NMAX, 2 )**2 > NWORK ) {
         INFO = -21;
      } else if ( 2*max( NMAX, 2 )**2 > LRWORK ) {
         INFO = -23;
      } else if ( 2*max( NMAX, 2 )**2 > LIWORK ) {
         INFO = -25;
      }

      if ( INFO != 0 ) {
         xerbla('ZDRVSG2STG', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) return;

      // More Important constants

      UNFL = dlamch( 'Safe minimum' );
      OVFL = dlamch( 'Overflow' );
      ULP = dlamch( 'Epsilon' )*dlamch( 'Base' );
      ULPINV = ONE / ULP;
      RTUNFL = sqrt( UNFL );
      RTOVFL = sqrt( OVFL );

      for (I = 1; I <= 4; I++) { // 20
         ISEED2[I] = ISEED( I );
      } // 20

      // Loop over sizes, types

      NERRS = 0;
      NMATS = 0;

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 650
         N = NN( JSIZE );
         ANINV = ONE / (max( 1, N )).toDouble();

         if ( NSIZES != 1 ) {
            MTYPES = min( MAXTYP, NTYPES );
         } else {
            MTYPES = min( MAXTYP+1, NTYPES );
         }

         KA9 = 0;
         KB9 = 0;
         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 640
            if( !DOTYPE( JTYPE ) ) GO TO 640;
            NMATS = NMATS + 1;
            NTEST = 0;

            for (J = 1; J <= 4; J++) { // 30
               IOLDSD[J] = ISEED( J );
            } // 30

            // 2)      Compute "A"

                    // Control parameters:

                // KMAGN  KMODE        KTYPE
            // =1  O(1)   clustered 1  zero
            // =2  large  clustered 2  identity
            // =3  small  exponential  (none)
            // =4         arithmetic   diagonal, w/ eigenvalues
            // =5         random log   hermitian, w/ eigenvalues
            // =6         random       (none)
            // =7                      random diagonal
            // =8                      random hermitian
            // =9                      banded, w/ eigenvalues

            if (MTYPES > MAXTYP) GO TO 90;

            ITYPE = KTYPE( JTYPE );
            IMODE = KMODE( JTYPE );

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE );

            } // 40
            ANORM = ONE;
            GO TO 70;

            } // 50
            ANORM = ( RTOVFL*ULP )*ANINV;
            GO TO 70;

            } // 60
            ANORM = RTUNFL*N*ULPINV;
            GO TO 70;

            } // 70

            IINFO = 0;
            COND = ULPINV;

            // Special Matrices -- Identity & Jordan block

            if ( ITYPE == 1 ) {

               // Zero

               KA = 0;
               KB = 0;
               zlaset('Full', LDA, N, CZERO, CZERO, A, LDA );

            } else if ( ITYPE == 2 ) {

               // Identity

               KA = 0;
               KB = 0;
               zlaset('Full', LDA, N, CZERO, CZERO, A, LDA );
               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A[JCOL, JCOL] = ANORM;
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               KA = 0;
               KB = 0;
               zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE == 5 ) {

               // Hermitian, eigenvalues specified

               KA = max( 0, N-1 );
               KB = KA;
               zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               KA = 0;
               KB = 0;
               zlatmr(N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // Hermitian, random eigenvalues

               KA = max( 0, N-1 );
               KB = KA;
               zlatmr(N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // Hermitian banded, eigenvalues specified

               // The following values are used for the half-bandwidths:

                 // ka = 1   kb = 1
                 // ka = 2   kb = 1
                 // ka = 2   kb = 2
                 // ka = 3   kb = 1
                 // ka = 3   kb = 2
                 // ka = 3   kb = 3

               KB9 = KB9 + 1;
               if ( KB9 > KA9 ) {
                  KA9 = KA9 + 1;
                  KB9 = 1;
               }
               KA = max( 0, min( N-1, KA9 ) );
               KB = max( 0, min( N-1, KB9 ) );
               zlatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, KA, KA, 'N', A, LDA, WORK, IINFO );

            } else {

               IINFO = 1;
            }

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD;
               INFO = ( IINFO ).abs();
               return;
            }

            } // 90

            ABSTOL = UNFL + UNFL;
            if ( N <= 1 ) {
               IL = 1;
               IU = N;
            } else {
               IL = 1 + INT( ( N-1 )*dlarnd( 1, ISEED2 ) );
               IU = 1 + INT( ( N-1 )*dlarnd( 1, ISEED2 ) );
               if ( IL > IU ) {
                  ITEMP = IL;
                  IL = IU;
                  IU = ITEMP;
               }
            }

            // 3) Call ZHEGV, ZHPGV, ZHBGV, CHEGVD, CHPGVD, CHBGVD,
               // ZHEGVX, ZHPGVX and ZHBGVX, do tests.

            // loop over the three generalized problems
                  // IBTYPE = 1: A*x = (lambda)*B*x
                  // IBTYPE = 2: A*B*x = (lambda)*x
                  // IBTYPE = 3: B*A*x = (lambda)*x

            for (IBTYPE = 1; IBTYPE <= 3; IBTYPE++) { // 630

               // loop over the setting UPLO

               for (IBUPLO = 1; IBUPLO <= 2; IBUPLO++) { // 620
                  if (IBUPLO == 1) UPLO = 'U';
                  IF( IBUPLO == 2 ) UPLO = 'L';

                  // Generate random well-conditioned positive definite
                  // matrix B, of bandwidth not greater than that of A.

                  zlatms(N, N, 'U', ISEED, 'P', RWORK, 5, TEN, ONE, KB, KB, UPLO, B, LDB, WORK( N+1 ), IINFO );

                  // Test ZHEGV

                  NTEST = NTEST + 1;

                  zlacpy(' ', N, N, A, LDA, Z, LDZ );
                  zlacpy(UPLO, N, N, B, LDB, BB, LDB );

                  zhegv(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, RWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHEGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                  zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  // Test ZHEGV_2STAGE

                  NTEST = NTEST + 1;

                  zlacpy(' ', N, N, A, LDA, Z, LDZ );
                  zlacpy(UPLO, N, N, B, LDB, BB, LDB );

                  zhegv_2stage(IBTYPE, 'N', UPLO, N, Z, LDZ, BB, LDB, D2, WORK, NWORK, RWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 ) 'ZHEGV_2STAGE(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                   // CALL ZSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z,
      // $                         LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  // Do Tests | D1 - D2 | / ( |D1| ulp )
                  // D1 computed using the standard 1-stage reduction as reference
                  // D2 computed using the 2-stage reduction

                  TEMP1 = ZERO;
                  TEMP2 = ZERO;
                  for (J = 1; J <= N; J++) { // 151
                     TEMP1 = max( TEMP1, ( D( J ) ).abs(),  ( D2( J ) ) ).abs();
                     TEMP2 = max( TEMP2, ABS( D( J )-D2( J ) ) );
                  } // 151

                  RESULT[NTEST] = TEMP2 /  max( UNFL, ULP*max( TEMP1, TEMP2 ) );

                  // Test ZHEGVD

                  NTEST = NTEST + 1;

                  zlacpy(' ', N, N, A, LDA, Z, LDZ );
                  zlacpy(UPLO, N, N, B, LDB, BB, LDB );

                  zhegvd(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHEGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                  zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  // Test ZHEGVX

                  NTEST = NTEST + 1;

                  zlacpy(' ', N, N, A, LDA, AB, LDA );
                  zlacpy(UPLO, N, N, B, LDB, BB, LDB );

                  zhegvx(IBTYPE, 'V', 'A', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHEGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                  zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1;

                  zlacpy(' ', N, N, A, LDA, AB, LDA );
                  zlacpy(UPLO, N, N, B, LDB, BB, LDB );

                  // since we do not know the exact eigenvalues of this
                  // eigenpair, we just set VL and VU as constants.
                  // It is quite possible that there are no eigenvalues
                  // in this interval.

                  VL = ZERO;
                  VU = ANORM;
                  zhegvx(IBTYPE, 'V', 'V', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHEGVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                  zsgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1;

                  zlacpy(' ', N, N, A, LDA, AB, LDA );
                  zlacpy(UPLO, N, N, B, LDB, BB, LDB );

                  zhegvx(IBTYPE, 'V', 'I', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHEGVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                  zsgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  } // 100

                  // Test ZHPGV

                  NTEST = NTEST + 1;

                  // Copy the matrices into packed storage.

                  if ( lsame( UPLO, 'U' ) ) {
                     IJ = 1;
                     for (J = 1; J <= N; J++) { // 120
                        for (I = 1; I <= J; I++) { // 110
                           AP[IJ] = A( I, J );
                           BP[IJ] = B( I, J );
                           IJ = IJ + 1;
                        } // 110
                     } // 120
                  } else {
                     IJ = 1;
                     for (J = 1; J <= N; J++) { // 140
                        for (I = J; I <= N; I++) { // 130
                           AP[IJ] = A( I, J );
                           BP[IJ] = B( I, J );
                           IJ = IJ + 1;
                        } // 130
                     } // 140
                  }

                  zhpgv(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, RWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHPGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 310;
                     }
                  }

                  // Do Test

                  zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  // Test ZHPGVD

                  NTEST = NTEST + 1;

                  // Copy the matrices into packed storage.

                  if ( lsame( UPLO, 'U' ) ) {
                     IJ = 1;
                     for (J = 1; J <= N; J++) { // 160
                        for (I = 1; I <= J; I++) { // 150
                           AP[IJ] = A( I, J );
                           BP[IJ] = B( I, J );
                           IJ = IJ + 1;
                        } // 150
                     } // 160
                  } else {
                     IJ = 1;
                     for (J = 1; J <= N; J++) { // 180
                        for (I = J; I <= N; I++) { // 170
                           AP[IJ] = A( I, J );
                           BP[IJ] = B( I, J );
                           IJ = IJ + 1;
                        } // 170
                     } // 180
                  }

                  zhpgvd(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHPGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 310;
                     }
                  }

                  // Do Test

                  zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  // Test ZHPGVX

                  NTEST = NTEST + 1;

                  // Copy the matrices into packed storage.

                  if ( lsame( UPLO, 'U' ) ) {
                     IJ = 1;
                     for (J = 1; J <= N; J++) { // 200
                        for (I = 1; I <= J; I++) { // 190
                           AP[IJ] = A( I, J );
                           BP[IJ] = B( I, J );
                           IJ = IJ + 1;
                        } // 190
                     } // 200
                  } else {
                     IJ = 1;
                     for (J = 1; J <= N; J++) { // 220
                        for (I = J; I <= N; I++) { // 210
                           AP[IJ] = A( I, J );
                           BP[IJ] = B( I, J );
                           IJ = IJ + 1;
                        } // 210
                     } // 220
                  }

                  zhpgvx(IBTYPE, 'V', 'A', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHPGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 310;
                     }
                  }

                  // Do Test

                  zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1;

                  // Copy the matrices into packed storage.

                  if ( lsame( UPLO, 'U' ) ) {
                     IJ = 1;
                     for (J = 1; J <= N; J++) { // 240
                        for (I = 1; I <= J; I++) { // 230
                           AP[IJ] = A( I, J );
                           BP[IJ] = B( I, J );
                           IJ = IJ + 1;
                        } // 230
                     } // 240
                  } else {
                     IJ = 1;
                     for (J = 1; J <= N; J++) { // 260
                        for (I = J; I <= N; I++) { // 250
                           AP[IJ] = A( I, J );
                           BP[IJ] = B( I, J );
                           IJ = IJ + 1;
                        } // 250
                     } // 260
                  }

                  VL = ZERO;
                  VU = ANORM;
                  zhpgvx(IBTYPE, 'V', 'V', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHPGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 310;
                     }
                  }

                  // Do Test

                  zsgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1;

                  // Copy the matrices into packed storage.

                  if ( lsame( UPLO, 'U' ) ) {
                     IJ = 1;
                     for (J = 1; J <= N; J++) { // 280
                        for (I = 1; I <= J; I++) { // 270
                           AP[IJ] = A( I, J );
                           BP[IJ] = B( I, J );
                           IJ = IJ + 1;
                        } // 270
                     } // 280
                  } else {
                     IJ = 1;
                     for (J = 1; J <= N; J++) { // 300
                        for (I = J; I <= N; I++) { // 290
                           AP[IJ] = A( I, J );
                           BP[IJ] = B( I, J );
                           IJ = IJ + 1;
                        } // 290
                     } // 300
                  }

                  zhpgvx(IBTYPE, 'V', 'I', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'ZHPGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 310;
                     }
                  }

                  // Do Test

                  zsgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  } // 310

                  if ( IBTYPE == 1 ) {

                     // TEST ZHBGV

                     NTEST = NTEST + 1;

                     // Copy the matrices into band storage.

                     if ( lsame( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 340
                           for (I = max( 1, J-KA ); I <= J; I++) { // 320
                              AB[KA+1+I-J, J] = A( I, J );
                           } // 320
                           for (I = max( 1, J-KB ); I <= J; I++) { // 330
                              BB[KB+1+I-J, J] = B( I, J );
                           } // 330
                        } // 340
                     } else {
                        for (J = 1; J <= N; J++) { // 370
                           for (I = J; I <= min( N, J+KA ); I++) { // 350
                              AB[1+I-J, J] = A( I, J );
                           } // 350
                           for (I = J; I <= min( N, J+KB ); I++) { // 360
                              BB[1+I-J, J] = B( I, J );
                           } // 360
                        } // 370
                     }

                     zhbgv('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, RWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'ZHBGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[NTEST] = ULPINV;
                           GO TO 620;
                        }
                     }

                     // Do Test

                     zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                     // TEST ZHBGVD

                     NTEST = NTEST + 1;

                     // Copy the matrices into band storage.

                     if ( lsame( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 400
                           for (I = max( 1, J-KA ); I <= J; I++) { // 380
                              AB[KA+1+I-J, J] = A( I, J );
                           } // 380
                           for (I = max( 1, J-KB ); I <= J; I++) { // 390
                              BB[KB+1+I-J, J] = B( I, J );
                           } // 390
                        } // 400
                     } else {
                        for (J = 1; J <= N; J++) { // 430
                           for (I = J; I <= min( N, J+KA ); I++) { // 410
                              AB[1+I-J, J] = A( I, J );
                           } // 410
                           for (I = J; I <= min( N, J+KB ); I++) { // 420
                              BB[1+I-J, J] = B( I, J );
                           } // 420
                        } // 430
                     }

                     zhbgvd('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'ZHBGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[NTEST] = ULPINV;
                           GO TO 620;
                        }
                     }

                     // Do Test

                     zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                     // Test ZHBGVX

                     NTEST = NTEST + 1;

                     // Copy the matrices into band storage.

                     if ( lsame( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 460
                           for (I = max( 1, J-KA ); I <= J; I++) { // 440
                              AB[KA+1+I-J, J] = A( I, J );
                           } // 440
                           for (I = max( 1, J-KB ); I <= J; I++) { // 450
                              BB[KB+1+I-J, J] = B( I, J );
                           } // 450
                        } // 460
                     } else {
                        for (J = 1; J <= N; J++) { // 490
                           for (I = J; I <= min( N, J+KA ); I++) { // 470
                              AB[1+I-J, J] = A( I, J );
                           } // 470
                           for (I = J; I <= min( N, J+KB ); I++) { // 480
                              BB[1+I-J, J] = B( I, J );
                           } // 480
                        } // 490
                     }

                     zhbgvx('V', 'A', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, max( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'ZHBGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[NTEST] = ULPINV;
                           GO TO 620;
                        }
                     }

                     // Do Test

                     zsgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                     NTEST = NTEST + 1;

                     // Copy the matrices into band storage.

                     if ( lsame( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 520
                           for (I = max( 1, J-KA ); I <= J; I++) { // 500
                              AB[KA+1+I-J, J] = A( I, J );
                           } // 500
                           for (I = max( 1, J-KB ); I <= J; I++) { // 510
                              BB[KB+1+I-J, J] = B( I, J );
                           } // 510
                        } // 520
                     } else {
                        for (J = 1; J <= N; J++) { // 550
                           for (I = J; I <= min( N, J+KA ); I++) { // 530
                              AB[1+I-J, J] = A( I, J );
                           } // 530
                           for (I = J; I <= min( N, J+KB ); I++) { // 540
                              BB[1+I-J, J] = B( I, J );
                           } // 540
                        } // 550
                     }

                     VL = ZERO;
                     VU = ANORM;
                     zhbgvx('V', 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, max( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'ZHBGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[NTEST] = ULPINV;
                           GO TO 620;
                        }
                     }

                     // Do Test

                     zsgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                     NTEST = NTEST + 1;

                     // Copy the matrices into band storage.

                     if ( lsame( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 580
                           for (I = max( 1, J-KA ); I <= J; I++) { // 560
                              AB[KA+1+I-J, J] = A( I, J );
                           } // 560
                           for (I = max( 1, J-KB ); I <= J; I++) { // 570
                              BB[KB+1+I-J, J] = B( I, J );
                           } // 570
                        } // 580
                     } else {
                        for (J = 1; J <= N; J++) { // 610
                           for (I = J; I <= min( N, J+KA ); I++) { // 590
                              AB[1+I-J, J] = A( I, J );
                           } // 590
                           for (I = J; I <= min( N, J+KB ); I++) { // 600
                              BB[1+I-J, J] = B( I, J );
                           } // 600
                        } // 610
                     }

                     zhbgvx('V', 'I', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, max( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'ZHBGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[NTEST] = ULPINV;
                           GO TO 620;
                        }
                     }

                     // Do Test

                     zsgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  }

               } // 620
            } // 630

            // End of Loop -- Check for RESULT(j) > THRESH

            NTESTT = NTESTT + NTEST;
            dlafts('ZSG', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS );
         } // 640
      } // 650

      // Summary

      dlasum('ZSG', NOUNIT, NERRS, NTESTT );

      return;

 9999 FORMAT( ' ZDRVSG2STG: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' );
      }
