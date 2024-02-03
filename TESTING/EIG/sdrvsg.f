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
      const              ZERO = 0.0E0, ONE = 1.0E0, TEN = 10.0E0 ;
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
      DATA               KMAGN / 2*1, 1, 1, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 6*1 /       DATA               KMODE / 2*0, 4, 3, 1, 4, 4, 4, 3, 1, 4, 4, 0, 0, 0, 6*4 /
      // ..
      // .. Executable Statements ..

      // 1)      Check for errors

      NTESTT = 0
      INFO = 0

      BADNN = .FALSE.
      NMAX = 0
      DO 10 J = 1, NSIZES
         NMAX = MAX( NMAX, NN( J ) )
         IF( NN( J ).LT.0 ) BADNN = .TRUE.
   10 CONTINUE

      // Check for errors

      if ( NSIZES.LT.0 ) {
         INFO = -1
      } else if ( BADNN ) {
         INFO = -2
      } else if ( NTYPES.LT.0 ) {
         INFO = -3
      } else if ( LDA.LE.1 .OR. LDA.LT.NMAX ) {
         INFO = -9
      } else if ( LDZ.LE.1 .OR. LDZ.LT.NMAX ) {
         INFO = -16
      } else if ( 2*MAX( NMAX, 3 )**2.GT.NWORK ) {
         INFO = -21
      } else if ( 2*MAX( NMAX, 3 )**2.GT.LIWORK ) {
         INFO = -23
      }

      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SDRVSG', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( NSIZES.EQ.0 .OR. NTYPES.EQ.0 ) RETURN

      // More Important constants

      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = SLAMCH( 'Overflow' )
      ULP = SLAMCH( 'Epsilon' )*SLAMCH( 'Base' )
      ULPINV = ONE / ULP
      RTUNFL = SQRT( UNFL )
      RTOVFL = SQRT( OVFL )

      DO 20 I = 1, 4
         ISEED2( I ) = ISEED( I )
   20 CONTINUE

      // Loop over sizes, types

      NERRS = 0
      NMATS = 0

      DO 650 JSIZE = 1, NSIZES
         N = NN( JSIZE )
         ANINV = ONE / REAL( MAX( 1, N ) )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         KA9 = 0
         KB9 = 0
         DO 640 JTYPE = 1, MTYPES
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 640
            NMATS = NMATS + 1
            NTEST = 0

            DO 30 J = 1, 4
               IOLDSD( J ) = ISEED( J )
   30       CONTINUE

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

            IF( MTYPES.GT.MAXTYP ) GO TO 90

            ITYPE = KTYPE( JTYPE )
            IMODE = KMODE( JTYPE )

            // Compute norm

            GO TO ( 40, 50, 60 )KMAGN( JTYPE )

   40       CONTINUE
            ANORM = ONE
            GO TO 70

   50       CONTINUE
            ANORM = ( RTOVFL*ULP )*ANINV
            GO TO 70

   60       CONTINUE
            ANORM = RTUNFL*N*ULPINV
            GO TO 70

   70       CONTINUE

            IINFO = 0
            COND = ULPINV

            // Special Matrices -- Identity & Jordan block

            if ( ITYPE.EQ.1 ) {

               // Zero

               KA = 0
               KB = 0
               CALL SLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )

            } else if ( ITYPE.EQ.2 ) {

               // Identity

               KA = 0
               KB = 0
               CALL SLASET( 'Full', LDA, N, ZERO, ZERO, A, LDA )
               DO 80 JCOL = 1, N
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE

            } else if ( ITYPE.EQ.4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               KA = 0
               KB = 0
               CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK( N+1 ), IINFO )

            } else if ( ITYPE.EQ.5 ) {

               // symmetric, eigenvalues specified

               KA = MAX( 0, N-1 )
               KB = KA
               CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK( N+1 ), IINFO )

            } else if ( ITYPE.EQ.7 ) {

               // Diagonal, random eigenvalues

               KA = 0
               KB = 0
               CALL SLATMR( N, N, 'S', ISEED, 'S', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.8 ) {

               // symmetric, random eigenvalues

               KA = MAX( 0, N-1 )
               KB = KA
               CALL SLATMR( N, N, 'S', ISEED, 'H', WORK, 6, ONE, ONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO )

            } else if ( ITYPE.EQ.9 ) {

               // symmetric banded, eigenvalues specified

               // The following values are used for the half-bandwidths:

                 // ka = 1   kb = 1
                 // ka = 2   kb = 1
                 // ka = 2   kb = 2
                 // ka = 3   kb = 1
                 // ka = 3   kb = 2
                 // ka = 3   kb = 3

               KB9 = KB9 + 1
               if ( KB9.GT.KA9 ) {
                  KA9 = KA9 + 1
                  KB9 = 1
               }
               KA = MAX( 0, MIN( N-1, KA9 ) )
               KB = MAX( 0, MIN( N-1, KB9 ) )
               CALL SLATMS( N, N, 'S', ISEED, 'S', WORK, IMODE, COND, ANORM, KA, KA, 'N', A, LDA, WORK( N+1 ), IINFO )

            } else {

               IINFO = 1
            }

            if ( IINFO.NE.0 ) {
               WRITE( NOUNIT, FMT = 9999 )'Generator', IINFO, N, JTYPE, IOLDSD
               INFO = ABS( IINFO )
               RETURN
            }

   90       CONTINUE

            ABSTOL = UNFL + UNFL
            if ( N.LE.1 ) {
               IL = 1
               IU = N
            } else {
               IL = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
               IU = 1 + INT( ( N-1 )*SLARND( 1, ISEED2 ) )
               if ( IL.GT.IU ) {
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

            DO 630 IBTYPE = 1, 3

               // loop over the setting UPLO

               DO 620 IBUPLO = 1, 2
                  IF( IBUPLO.EQ.1 ) UPLO = 'U'                   IF( IBUPLO.EQ.2 ) UPLO = 'L'

                  // Generate random well-conditioned positive definite
                  // matrix B, of bandwidth not greater than that of A.

                  CALL SLATMS( N, N, 'U', ISEED, 'P', WORK, 5, TEN, ONE, KB, KB, UPLO, B, LDB, WORK( N+1 ), IINFO )

                  // Test SSYGV

                  NTEST = NTEST + 1

                  CALL SLACPY( ' ', N, N, A, LDA, Z, LDZ )
                  CALL SLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  CALL SSYGV( IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  CALL SSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                  // Test SSYGVD

                  NTEST = NTEST + 1

                  CALL SLACPY( ' ', N, N, A, LDA, Z, LDZ )
                  CALL SLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  CALL SSYGVD( IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, IWORK, LIWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  CALL SSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                  // Test SSYGVX

                  NTEST = NTEST + 1

                  CALL SLACPY( ' ', N, N, A, LDA, AB, LDA )
                  CALL SLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  CALL SSYGVX( IBTYPE, 'V', 'A', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK( N+1 ), IWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  CALL SSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                  NTEST = NTEST + 1

                  CALL SLACPY( ' ', N, N, A, LDA, AB, LDA )
                  CALL SLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  // since we do not know the exact eigenvalues of this
                  // eigenpair, we just set VL and VU as constants.
                  // It is quite possible that there are no eigenvalues
                  // in this interval.

                  VL = ZERO
                  VU = ANORM
                  CALL SSYGVX( IBTYPE, 'V', 'V', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK( N+1 ), IWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  CALL SSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                  NTEST = NTEST + 1

                  CALL SLACPY( ' ', N, N, A, LDA, AB, LDA )
                  CALL SLACPY( UPLO, N, N, B, LDB, BB, LDB )

                  CALL SSYGVX( IBTYPE, 'V', 'I', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, IWORK( N+1 ), IWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSYGVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  CALL SSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

  100             CONTINUE

                  // Test SSPGV

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     DO 120 J = 1, N
                        DO 110 I = 1, J
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  110                   CONTINUE
  120                CONTINUE
                  } else {
                     IJ = 1
                     DO 140 J = 1, N
                        DO 130 I = J, N
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  130                   CONTINUE
  140                CONTINUE
                  }

                  CALL SSPGV( IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  CALL SSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                  // Test SSPGVD

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     DO 160 J = 1, N
                        DO 150 I = 1, J
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  150                   CONTINUE
  160                CONTINUE
                  } else {
                     IJ = 1
                     DO 180 J = 1, N
                        DO 170 I = J, N
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  170                   CONTINUE
  180                CONTINUE
                  }

                  CALL SSPGVD( IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, NWORK, IWORK, LIWORK, IINFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  CALL SSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                  // Test SSPGVX

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     DO 200 J = 1, N
                        DO 190 I = 1, J
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  190                   CONTINUE
  200                CONTINUE
                  } else {
                     IJ = 1
                     DO 220 J = 1, N
                        DO 210 I = J, N
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  210                   CONTINUE
  220                CONTINUE
                  }

                  CALL SSPGVX( IBTYPE, 'V', 'A', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, INFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  CALL SSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     DO 240 J = 1, N
                        DO 230 I = 1, J
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  230                   CONTINUE
  240                CONTINUE
                  } else {
                     IJ = 1
                     DO 260 J = 1, N
                        DO 250 I = J, N
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  250                   CONTINUE
  260                CONTINUE
                  }

                  VL = ZERO
                  VU = ANORM
                  CALL SSPGVX( IBTYPE, 'V', 'V', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, INFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  CALL SSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     DO 280 J = 1, N
                        DO 270 I = 1, J
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  270                   CONTINUE
  280                CONTINUE
                  } else {
                     IJ = 1
                     DO 300 J = 1, N
                        DO 290 I = J, N
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  290                   CONTINUE
  300                CONTINUE
                  }

                  CALL SSPGVX( IBTYPE, 'V', 'I', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, INFO )
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'SSPGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  CALL SSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

  310             CONTINUE

                  if ( IBTYPE.EQ.1 ) {

                     // TEST SSBGV

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        DO 340 J = 1, N
                           DO 320 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  320                      CONTINUE
                           DO 330 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  330                      CONTINUE
  340                   CONTINUE
                     } else {
                        DO 370 J = 1, N
                           DO 350 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  350                      CONTINUE
                           DO 360 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  360                      CONTINUE
  370                   CONTINUE
                     }

                     CALL SSBGV( 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, IINFO )
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     CALL SSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                     // TEST SSBGVD

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        DO 400 J = 1, N
                           DO 380 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  380                      CONTINUE
                           DO 390 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  390                      CONTINUE
  400                   CONTINUE
                     } else {
                        DO 430 J = 1, N
                           DO 410 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  410                      CONTINUE
                           DO 420 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  420                      CONTINUE
  430                   CONTINUE
                     }

                     CALL SSBGVD( 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, NWORK, IWORK, LIWORK, IINFO )
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     CALL SSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                     // Test SSBGVX

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        DO 460 J = 1, N
                           DO 440 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  440                      CONTINUE
                           DO 450 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  450                      CONTINUE
  460                   CONTINUE
                     } else {
                        DO 490 J = 1, N
                           DO 470 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  470                      CONTINUE
                           DO 480 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  480                      CONTINUE
  490                   CONTINUE
                     }

                     CALL SSBGVX( 'V', 'A', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, IINFO )
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     CALL SSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )


                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        DO 520 J = 1, N
                           DO 500 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  500                      CONTINUE
                           DO 510 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  510                      CONTINUE
  520                   CONTINUE
                     } else {
                        DO 550 J = 1, N
                           DO 530 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  530                      CONTINUE
                           DO 540 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  540                      CONTINUE
  550                   CONTINUE
                     }

                     VL = ZERO
                     VU = ANORM
                     CALL SSBGVX( 'V', 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, IINFO )
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     CALL SSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        DO 580 J = 1, N
                           DO 560 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  560                      CONTINUE
                           DO 570 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  570                      CONTINUE
  580                   CONTINUE
                     } else {
                        DO 610 J = 1, N
                           DO 590 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  590                      CONTINUE
                           DO 600 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  600                      CONTINUE
  610                   CONTINUE
                     }

                     CALL SSBGVX( 'V', 'I', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, IWORK( N+1 ), IWORK, IINFO )
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'SSBGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     CALL SSGT01( IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RESULT( NTEST ) )

                  }

  620          CONTINUE
  630       CONTINUE

            // End of Loop -- Check for RESULT(j) > THRESH

            NTESTT = NTESTT + NTEST
            CALL SLAFTS( 'SSG', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS )
  640    CONTINUE
  650 CONTINUE

      // Summary

      CALL SLASUM( 'SSG', NOUNIT, NERRS, NTESTT )

      RETURN

      // End of SDRVSG

 9999 FORMAT( ' SDRVSG: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )
      }
