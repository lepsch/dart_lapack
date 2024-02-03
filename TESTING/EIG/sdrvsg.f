      SUBROUTINE SDRVSG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, LDB, D, Z, LDZ, AB, BB, AP, BP, WORK, NWORK, IWORK, LIWORK, RESULT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LDZ, LIWORK, NOUNIT, NSIZES, NTYPES, NWORK;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      REAL               A( LDA, * ), AB( LDA, * ), AP( * ), B( LDB, * ), BB( LDB, * ), BP( * ), D( * ), RESULT( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TEN
      const              ZERO = 0.0, ONE = 1.0, TEN = 10.0 ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             UPLO;
      int                I, IBTYPE, IBUPLO, IINFO, IJ, IL, IMODE, ITEMP, ITYPE, IU, J, JCOL, JSIZE, JTYPE, KA, KA9, KB, KB9, M, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      REAL               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU
      // ..
      // .. Local Arrays ..
      int                IDUMMA( 1 ), IOLDSD( 4 ), ISEED2( 4 ), KMAGN( MAXTYP ), KMODE( MAXTYP ), KTYPE( MAXTYP );
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, SLARND
      // EXTERNAL LSAME, SLAMCH, SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACPY, SLAFTS, SLASET, SLASUM, SLATMR, SLATMS, SSBGV, SSBGVD, SSBGVX, SSGT01, SSPGV, SSPGVD, SSPGVX, SSYGV, SSYGVD, SSYGVX, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, REAL, SQRT
      // ..
      // .. Data statements ..
      DATA               KTYPE / 1, 2, 5*4, 5*5, 3*8, 6*9 /
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 6*1 /
      DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 6*4 /
      // ..
      // .. Executable Statements ..

      // 1)      Check for errors

      NTESTT = 0
      INFO = 0

      BADNN = false;
      NMAX = 0
      for (J = 1; J <= NSIZES; J++) { // 10
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ) < 0 ) BADNN = true;
      } // 10

      // Check for errors

      if ( NSIZES < 0 ) {
         INFO = -1
      } else if ( BADNN ) {
         INFO = -2
      } else if ( NTYPES < 0 ) {
         INFO = -3
      } else if ( LDA <= 1 || LDA < NMAX ) {
         INFO = -9
      } else if ( LDZ <= 1 || LDZ < NMAX ) {
         INFO = -16
      } else if ( 2*MAX( NMAX, 3 )**2 > NWORK ) {
         INFO = -21
      } else if ( 2*MAX( NMAX, 3 )**2 > LIWORK ) {
         INFO = -23
      }

      if ( INFO != 0 ) {
         xerbla('SDRVSG', -INFO );
         RETURN
      }

      // Quick return if possible

      if (NSIZES == 0 || NTYPES == 0) RETURN;

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = SLAMCH( 'Overflow' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )

      for (I = 1; I <= 4; I++) { // 20
         ISEED2( I ) = ISEED( I )
      } // 20

      // Loop over sizes, types

      NERRS = 0
      NMATS = 0

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 650
         N = NN( JSIZE )
         ANINV = ONE / REAL( MAX( 1, N ) )

         if ( NSIZES != 1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         KA9 = 0
         KB9 = 0
         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 640
            IF( !DOTYPE( JTYPE ) ) GO TO 640
            NMATS = NMATS + 1
            NTEST = 0

            for (J = 1; J <= 4; J++) { // 30
               IOLDSD( J ) = ISEED( J )
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

            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE )

            } // 40
            ANORM = ONE
            GO TO 70

            } // 50
            ANORM = ( RTOVFL*ULP )*ANINV
            GO TO 70

            } // 60
            ANORM = RTUNFL*N*ULPINV
            GO TO 70

            } // 70

            IINFO = 0
            COND = ULPINV

            // Special Matrices -- Identity & Jordan block

            if ( ITYPE == 1 ) {

               // Zero

               KA = 0
               KB = 0
               slaset('Full', LDA, N, ZERO, ZERO, A, LDA );

            } else if ( ITYPE == 2 ) {

               // Identity

               KA = 0
               KB = 0
               slaset('Full', LDA, N, ZERO, ZERO, A, LDA );
               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A( JCOL, JCOL ) = ANORM
               } // 80

            } else if ( ITYPE == 4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               KA = 0
               KB = 0
               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 5 ) {

               // symmetric, eigenvalues specified

               KA = MAX( 0, N-1 )
               KB = KA
               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else if ( ITYPE == 7 ) {

               // Diagonal, random eigenvalues

               KA = 0
               KB = 0
               slatmr(N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE == 8 ) {

               // symmetric, random eigenvalues

               KA = MAX( 0, N-1 )
               KB = KA
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

               KB9 = KB9 + 1
               if ( KB9 > KA9 ) {
                  KA9 = KA9 + 1
                  KB9 = 1
               }
               KA = MAX( 0, MIN( N-1, KA9 ) )
               KB = MAX( 0, MIN( N-1, KB9 ) )
               slatms(N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, KA, KA, 'N', A, LDA, WORK( N+1 ), IINFO );

            } else {

               IINFO = 1
            }

            if ( IINFO != 0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

            } // 90

            ABSTOL = UNFL + UNFL
            if ( N <= 1 ) {
               IL = 1
               IU = N
            } else {
               IL = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
               IU = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
               if ( IL > IU ) {
                  ITEMP = IL
                  IL = IU
                  IU = ITEMP
               }
            }

            // 3) Call SSYGV, SSPGV, SSBGV, SSYGVD, SSPGVD, SSBGVD,
               // SSYGVX, SSPGVX, and SSBGVX, do tests.

            // loop over the three generalized problems
                  // IBTYPE = 1: A*x = (lambda)*B*x
                  // IBTYPE = 2: A*B*x = (lambda)*x
                  // IBTYPE = 3: B*A*x = (lambda)*x

            for (IBTYPE = 1; IBTYPE <= 3; IBTYPE++) { // 630

               // loop over the setting UPLO

               for (IBUPLO = 1; IBUPLO <= 2; IBUPLO++) { // 620
                  if (IBUPLO == 1) UPLO = 'U'                   IF( IBUPLO == 2 ) UPLO = 'L';

                  // Generate random well-conditioned positive definite
                  // matrix B, of bandwidth not greater than that of A.

                  slatms(N, N, 'U', ISEED, 'P', WORK, 5, TEN, ONE, KB, KB, UPLO, B, LDB, WORK( N+1 ), IINFO );

                  // Test SSYGV

                  NTEST = NTEST + 1

                  slacpy(' ', N, N, A, LDA, Z, LDZ );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  ssygv(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO < 0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  // Test SSYGVD

                  NTEST = NTEST + 1

                  slacpy(' ', N, N, A, LDA, Z, LDZ );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  ssygvd(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, IWORK, LIWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO < 0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  // Test SSYGVX

                  NTEST = NTEST + 1

                  slacpy(' ', N, N, A, LDA, AB, LDA );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  ssygvx(IBTYPE, 'V', 'A', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO < 0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1

                  slacpy(' ', N, N, A, LDA, AB, LDA );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  // since we do not know the exact eigenvalues of this
                  // eigenpair, we just set VL and VU as constants.
                  // It is quite possible that there are no eigenvalues
                  // in this interval.

                  VL = ZERO
                  VU = ANORM
                  ssygvx(IBTYPE, 'V', 'V', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO < 0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1

                  slacpy(' ', N, N, A, LDA, AB, LDA );
                  slacpy(UPLO, N, N, B, LDB, BB, LDB );

                  ssygvx(IBTYPE, 'V', 'I', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO < 0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  } // 100

                  // Test SSPGV

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 120
                        for (I = 1; I <= J; I++) { // 110
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
                        } // 110
                     } // 120
                  } else {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 140
                        for (I = J; I <= N; I++) { // 130
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
                        } // 130
                     } // 140
                  }

                  sspgv(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO < 0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  // Test SSPGVD

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 160
                        for (I = 1; I <= J; I++) { // 150
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
                        } // 150
                     } // 160
                  } else {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 180
                        for (I = J; I <= N; I++) { // 170
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
                        } // 170
                     } // 180
                  }

                  sspgvd(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, NWORK, IWORK, LIWORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO < 0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  // Test SSPGVX

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 200
                        for (I = 1; I <= J; I++) { // 190
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
                        } // 190
                     } // 200
                  } else {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 220
                        for (I = J; I <= N; I++) { // 210
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
                        } // 210
                     } // 220
                  }

                  sspgvx(IBTYPE, 'V', 'A', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO < 0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 240
                        for (I = 1; I <= J; I++) { // 230
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
                        } // 230
                     } // 240
                  } else {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 260
                        for (I = J; I <= N; I++) { // 250
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
                        } // 250
                     } // 260
                  }

                  VL = ZERO
                  VU = ANORM
                  sspgvx(IBTYPE, 'V', 'V', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO < 0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 280
                        for (I = 1; I <= J; I++) { // 270
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
                        } // 270
                     } // 280
                  } else {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 300
                        for (I = J; I <= N; I++) { // 290
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
                        } // 290
                     } // 300
                  }

                  sspgvx(IBTYPE, 'V', 'I', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO < 0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  } // 310

                  if ( IBTYPE == 1 ) {

                     // TEST SSBGV

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 340
                           DO 320 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
                           } // 320
                           DO 330 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
                           } // 330
                        } // 340
                     } else {
                        for (J = 1; J <= N; J++) { // 370
                           DO 350 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
                           } // 350
                           DO 360 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
                           } // 360
                        } // 370
                     }

                     ssbgv('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO < 0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                     // TEST SSBGVD

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 400
                           DO 380 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
                           } // 380
                           DO 390 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
                           } // 390
                        } // 400
                     } else {
                        for (J = 1; J <= N; J++) { // 430
                           DO 410 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
                           } // 410
                           DO 420 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
                           } // 420
                        } // 430
                     }

                     ssbgvd('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, NWORK, IWORK, LIWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO < 0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     ssgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                     // Test SSBGVX

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 460
                           DO 440 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
                           } // 440
                           DO 450 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
                           } // 450
                        } // 460
                     } else {
                        for (J = 1; J <= N; J++) { // 490
                           DO 470 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
                           } // 470
                           DO 480 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
                           } // 480
                        } // 490
                     }

                     ssbgvx('V', 'A', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO < 0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );


                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 520
                           DO 500 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
                           } // 500
                           DO 510 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
                           } // 510
                        } // 520
                     } else {
                        for (J = 1; J <= N; J++) { // 550
                           DO 530 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
                           } // 530
                           DO 540 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
                           } // 540
                        } // 550
                     }

                     VL = ZERO
                     VU = ANORM
                     ssbgvx('V', 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO < 0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 580
                           DO 560 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
                           } // 560
                           DO 570 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
                           } // 570
                        } // 580
                     } else {
                        for (J = 1; J <= N; J++) { // 610
                           DO 590 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
                           } // 590
                           DO 600 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
                           } // 600
                        } // 610
                     }

                     ssbgvx('V', 'I', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO != 0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO < 0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     ssgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) );

                  }

               } // 620
            } // 630

            // End of Loop -- Check for RESULT(j) > THRESH

            NTESTT = NTESTT + NTEST
            slafts('SSG', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS );
         } // 640
      } // 650

      // Summary

      slasum('SSG', NOUNIT, NERRS, NTESTT );

      RETURN

      // End of SDRVSG

 9999 FORMAT( ' SDRVSG: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
      }
