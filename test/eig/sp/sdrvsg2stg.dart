      void sdrvsg2stg(final int NSIZES, final int NN, final int NTYPES, final Array<bool> DOTYPE_, final Array<int> ISEED_, final int THRESH, final int NOUNIT, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int D, final int D2, final Matrix<double> Z_, final int LDZ, final int AB, final int BB, final int AP, final int BP, final Array<double> _WORK_, final int NWORK, final Array<int> IWORK_, final int LIWORK, final int RESULT, final Box<int> INFO,) {
  final DOTYPE = DOTYPE_.dim();
  final ISEED = ISEED_.dim();
  final A = A_.dim();
  final B = B_.dim();
  final Z = Z_.dim();
  final _WORK = _WORK_.dim();
  final IWORK = IWORK_.dim();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LDZ, LIWORK, NOUNIT, NSIZES, NTYPES, NWORK;
      double               THRESH;
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      double               A( LDA, * ), AB( LDA, * ), AP( * ), B( LDB, * ), BB( LDB, * ), BP( * ), D( * ), D2( * ), RESULT( * ), WORK( * ), Z( LDZ, * );
      // ..

      double               ZERO, ONE, TEN;
      const              ZERO = 0.0, ONE = 1.0, TEN = 10.0 ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      bool               BADNN;
      String             UPLO;
      int                I, IBTYPE, IBUPLO, IINFO, IJ, IL, IMODE, ITEMP, ITYPE, IU, J, JCOL, JSIZE, JTYPE, KA, KA9, KB, KB9, M, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      double               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU, TEMP1, TEMP2;
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               SLAMCH, SLARND;
      // EXTERNAL lsame, SLAMCH, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLAFTS, SLASET, SLASUM, SLATMR, SLATMS, SSBGV, SSBGVD, SSBGVX, SSGT01, SSPGV, SSPGVD, SSPGVX, SSYGV, SSYGVD, SSYGVX, XERBLA, SSYGV_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN, SQRT
      // ..
      // .. Data statements ..
      const KTYPE = [ 1, 2, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 8, 8, 8, 9, 9, 9, 9, 9, 9,];
      const KMAGN = [ 1, 1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 1, 1, 1,];
      const KMODE = [ 0, 0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 4, 4, 4, 4, 4, 4,];

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
      } else if ( 2*max( NMAX, 3 )**2 > NWORK ) {
         INFO = -21;
      } else if ( 2*max( NMAX, 3 )**2 > LIWORK ) {
         INFO = -23;
      }

      if ( INFO != 0 ) {
         xerbla('SDRVSG2STG', -INFO );
         return;
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) return;

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' );
      OVFL = SLAMCH( 'Overflow' );
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' );
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
         ANINV = ONE / REAL( max( 1, N ) );

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
               slaset('Full', LDA, N, ZERO, ZERO, A, LDA );

            } else if ( ITYPE == 2 ) {

               // Identity

               KA = 0;
               KB = 0;
               slaset('Full', LDA, N, ZERO, ZERO, A, LDA );
               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A[JCOL][JCOL] = ANORM;
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               KA = 0;
               KB = 0;
               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 5 ) {

               // symmetric, eigenvalues specified

               KA = max( 0, N-1 );
               KB = KA;
               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               KA = 0;
               KB = 0;
               slatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // symmetric, random eigenvalues

               KA = max( 0, N-1 );
               KB = KA;
               slatmr(N, N, 'S', ISEED, 'H', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 9 ) {

               // symmetric banded, eigenvalues specified

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
               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, KA, KA, 'N', A, LDA, WORK( N+1 ), IINFO );

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
               IL = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) );
               IU = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) );
               if ( IL > IU ) {
                  ITEMP = IL;
                  IL = IU;
                  IU = ITEMP;
               }
            }

            // 3) Call SSYGV, SSPGV, SSBGV, SSYGVD, SSPGVD, SSBGVD,
            //    SSYGVX, SSPGVX, and SSBGVX, do tests.

            // loop over the three generalized problems
            //       IBTYPE = 1: A*x = (lambda)*B*x
            //       IBTYPE = 2: A*B*x = (lambda)*x
            //       IBTYPE = 3: B*A*x = (lambda)*x

            for (IBTYPE = 1; IBTYPE <= 3; IBTYPE++) { // 630

               // loop over the setting UPLO

               for (IBUPLO = 1; IBUPLO <= 2; IBUPLO++) { // 620
                  if (IBUPLO == 1) UPLO = 'U';
                  IF( IBUPLO == 2 ) UPLO = 'L';

                  // Generate random well-conditioned positive definite
                  // matrix B, of bandwidth not greater than that of A.

                  slatms(N, N, 'U', ISEED, 'P', WORK, 5, TEN, ONE, KB, KB, UPLO, B, LDB, WORK( N+1 ), IINFO );

                  // Test SSYGV

                  NTEST = NTEST + 1;

                  slacpy(' ', N, N, A, LDA, Z, LDZ );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  ssygv(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  // Test SSYGV_2STAGE

                  NTEST = NTEST + 1;

                  slacpy(' ', N, N, A, LDA, Z, LDZ );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  ssygv_2stage(IBTYPE, 'N', UPLO, N, Z, LDZ, BB, LDB, D2, WORK, NWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 ) 'SSYGV_2STAGE(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                   // CALL SSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z,
      // $                         LDZ, D, WORK, RESULT( NTEST ) )


                  // Do Tests | D1 - D2 | / ( |D1| ulp )
                  // D1 computed using the standard 1-stage reduction as reference
                  // D2 computed using the 2-stage reduction

                  TEMP1 = ZERO;
                  TEMP2 = ZERO;
                  for (J = 1; J <= N; J++) { // 151
                     TEMP1 = max( TEMP1, ( D( J ) ).abs(),  ( D2( J ) ).abs() );
                     TEMP2 = max( TEMP2, ABS( D( J )-D2( J ) ) );
                  } // 151

                  RESULT[NTEST] = TEMP2 /  max( UNFL, ULP*max( TEMP1, TEMP2 ) );

                  // Test SSYGVD

                  NTEST = NTEST + 1;

                  slacpy(' ', N, N, A, LDA, Z, LDZ );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  ssygvd(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, IWORK, LIWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  // Test SSYGVX

                  NTEST = NTEST + 1;

                  slacpy(' ', N, N, A, LDA, AB, LDA );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  ssygvx(IBTYPE, 'V', 'A', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1;

                  slacpy(' ', N, N, A, LDA, AB, LDA );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  // since we do not know the exact eigenvalues of this
                  // eigenpair, we just set VL and VU as constants.
                  // It is quite possible that there are no eigenvalues
                  // in this interval.

                  VL = ZERO;
                  VU = ANORM;
                  ssygvx(IBTYPE, 'V', 'V', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1;

                  slacpy(' ', N, N, A, LDA, AB, LDA );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  ssygvx(IBTYPE, 'V', 'I', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 100;
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  } // 100

                  // Test SSPGV

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

                  sspgv(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 310;
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  // Test SSPGVD

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

                  sspgvd(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, NWORK, IWORK, LIWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 310;
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  // Test SSPGVX

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

                  sspgvx(IBTYPE, 'V', 'A', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 310;
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

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
                  sspgvx(IBTYPE, 'V', 'V', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 310;
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

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

                  sspgvx(IBTYPE, 'V', 'I', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                     INFO = ( IINFO ).abs();
                     if ( IINFO < 0 ) {
                        return;
                     } else {
                        RESULT[NTEST] = ULPINV;
                        GO TO 310;
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  } // 310

                  if ( IBTYPE == 1 ) {

                     // TEST SSBGV

                     NTEST = NTEST + 1;

                     // Copy the matrices into band storage.

                     if ( lsame( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 340
                           for (I = max( 1, J-KA ); I <= J; I++) { // 320
                              AB[KA+1+I-J][J] = A( I, J );
                           } // 320
                           for (I = max( 1, J-KB ); I <= J; I++) { // 330
                              BB[KB+1+I-J][J] = B( I, J );
                           } // 330
                        } // 340
                     } else {
                        for (J = 1; J <= N; J++) { // 370
                           for (I = J; I <= min( N, J+KA ); I++) { // 350
                              AB[1+I-J][J] = A( I, J );
                           } // 350
                           for (I = J; I <= min( N, J+KB ); I++) { // 360
                              BB[1+I-J][J] = B( I, J );
                           } // 360
                        } // 370
                     }

                     ssbgv('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[NTEST] = ULPINV;
                           GO TO 620;
                        }
                     }

                     // Do Test

                     ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                     // TEST SSBGVD

                     NTEST = NTEST + 1;

                     // Copy the matrices into band storage.

                     if ( lsame( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 400
                           for (I = max( 1, J-KA ); I <= J; I++) { // 380
                              AB[KA+1+I-J][J] = A( I, J );
                           } // 380
                           for (I = max( 1, J-KB ); I <= J; I++) { // 390
                              BB[KB+1+I-J][J] = B( I, J );
                           } // 390
                        } // 400
                     } else {
                        for (J = 1; J <= N; J++) { // 430
                           for (I = J; I <= min( N, J+KA ); I++) { // 410
                              AB[1+I-J][J] = A( I, J );
                           } // 410
                           for (I = J; I <= min( N, J+KB ); I++) { // 420
                              BB[1+I-J][J] = B( I, J );
                           } // 420
                        } // 430
                     }

                     ssbgvd('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, NWORK, IWORK, LIWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[NTEST] = ULPINV;
                           GO TO 620;
                        }
                     }

                     // Do Test

                     ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                     // Test SSBGVX

                     NTEST = NTEST + 1;

                     // Copy the matrices into band storage.

                     if ( lsame( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 460
                           for (I = max( 1, J-KA ); I <= J; I++) { // 440
                              AB[KA+1+I-J][J] = A( I, J );
                           } // 440
                           for (I = max( 1, J-KB ); I <= J; I++) { // 450
                              BB[KB+1+I-J][J] = B( I, J );
                           } // 450
                        } // 460
                     } else {
                        for (J = 1; J <= N; J++) { // 490
                           for (I = J; I <= min( N, J+KA ); I++) { // 470
                              AB[1+I-J][J] = A( I, J );
                           } // 470
                           for (I = J; I <= min( N, J+KB ); I++) { // 480
                              BB[1+I-J][J] = B( I, J );
                           } // 480
                        } // 490
                     }

                     ssbgvx('V', 'A', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, max( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[NTEST] = ULPINV;
                           GO TO 620;
                        }
                     }

                     // Do Test

                     ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );


                     NTEST = NTEST + 1;

                     // Copy the matrices into band storage.

                     if ( lsame( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 520
                           for (I = max( 1, J-KA ); I <= J; I++) { // 500
                              AB[KA+1+I-J][J] = A( I, J );
                           } // 500
                           for (I = max( 1, J-KB ); I <= J; I++) { // 510
                              BB[KB+1+I-J][J] = B( I, J );
                           } // 510
                        } // 520
                     } else {
                        for (J = 1; J <= N; J++) { // 550
                           for (I = J; I <= min( N, J+KA ); I++) { // 530
                              AB[1+I-J][J] = A( I, J );
                           } // 530
                           for (I = J; I <= min( N, J+KB ); I++) { // 540
                              BB[1+I-J][J] = B( I, J );
                           } // 540
                        } // 550
                     }

                     VL = ZERO;
                     VU = ANORM;
                     ssbgvx('V', 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, max( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[NTEST] = ULPINV;
                           GO TO 620;
                        }
                     }

                     // Do Test

                     ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                     NTEST = NTEST + 1;

                     // Copy the matrices into band storage.

                     if ( lsame( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 580
                           for (I = max( 1, J-KA ); I <= J; I++) { // 560
                              AB[KA+1+I-J][J] = A( I, J );
                           } // 560
                           for (I = max( 1, J-KB ); I <= J; I++) { // 570
                              BB[KB+1+I-J][J] = B( I, J );
                           } // 570
                        } // 580
                     } else {
                        for (J = 1; J <= N; J++) { // 610
                           for (I = J; I <= min( N, J+KA ); I++) { // 590
                              AB[1+I-J][J] = A( I, J );
                           } // 590
                           for (I = J; I <= min( N, J+KB ); I++) { // 600
                              BB[1+I-J][J] = B( I, J );
                           } // 600
                        } // 610
                     }

                     ssbgvx('V', 'I', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, max( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD;
                        INFO = ( IINFO ).abs();
                        if ( IINFO < 0 ) {
                           return;
                        } else {
                           RESULT[NTEST] = ULPINV;
                           GO TO 620;
                        }
                     }

                     // Do Test

                     ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  }

               } // 620
            } // 630

            // End of Loop -- Check for RESULT(j) > THRESH

            NTESTT = NTESTT + NTEST;
            slafts('SSG', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS );
         } // 640
      } // 650

      // Summary

      slasum('SSG', NOUNIT, NERRS, NTESTT );

      return;
 9999 FORMAT( ' SDRVSG2STG: ${} returned INFO=${.i6}.\n${' ' * 9}N=${.i6}, JTYPE=${.i6}, ISEED=(${.i5(4, ',')})' );
      }
