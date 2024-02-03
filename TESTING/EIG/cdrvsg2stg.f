      SUBROUTINE CDRVSG2STG( NSIZES, NN, NTYPES, DOTYPE, ISEED, THRESH, NOUNIT, A, LDA, B, LDB, D, D2, Z, LDZ, AB, BB, AP, BP, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, RESULT, INFO )

      IMPLICIT NONE

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LDZ, LIWORK, LRWORK, NOUNIT, NSIZES, NTYPES, NWORK;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                ISEED( 4 ), IWORK( * ), NN( * );
      REAL               D( * ), D2( * ), RESULT( * ), RWORK( * )
      COMPLEX            A( LDA, * ), AB( LDA, * ), AP( * ), B( LDB, * ), BB( LDB, * ), BP( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE, TEN
      const              ZERO = 0.0E+0, ONE = 1.0E+0, TEN = 10.0E+0 ;
      COMPLEX            CZERO, CONE
      const              CZERO = ( 0.0E+0, 0.0E+0 ), CONE = ( 1.0E+0, 0.0E+0 ) ;
      int                MAXTYP;
      const              MAXTYP = 21 ;
      // ..
      // .. Local Scalars ..
      bool               BADNN;
      String             UPLO;
      int                I, IBTYPE, IBUPLO, IINFO, IJ, IL, IMODE, ITEMP, ITYPE, IU, J, JCOL, JSIZE, JTYPE, KA, KA9, KB, KB9, M, MTYPES, N, NERRS, NMATS, NMAX, NTEST, NTESTT;
      REAL               ABSTOL, ANINV, ANORM, COND, OVFL, RTOVFL, RTUNFL, ULP, ULPINV, UNFL, VL, VU, TEMP1, TEMP2
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
      // EXTERNAL SLAFTS, SLASUM, XERBLA, CHBGV, CHBGVD, CHBGVX, CHEGV, CHEGVD, CHEGVX, CHPGV, CHPGVD, CHPGVX, CLACPY, CLASET, CLATMR, CLATMS, CSGT01, CHEGV_2STAGE
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, REAL, MAX, MIN, SQRT
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
      for (J = 1; J <= NSIZES; J++) { // 10
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
      } else if ( 2*MAX( NMAX, 2 )**2.GT.NWORK ) {
         INFO = -21
      } else if ( 2*MAX( NMAX, 2 )**2.GT.LRWORK ) {
         INFO = -23
      } else if ( 2*MAX( NMAX, 2 )**2.GT.LIWORK ) {
         INFO = -25
      }

      if ( INFO.NE.0 ) {
         xerbla('CDRVSG2STG', -INFO );
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

      for (I = 1; I <= 4; I++) { // 20
         ISEED2( I ) = ISEED( I )
   20 CONTINUE

      // Loop over sizes, types

      NERRS = 0
      NMATS = 0

      for (JSIZE = 1; JSIZE <= NSIZES; JSIZE++) { // 650
         N = NN( JSIZE )
         ANINV = ONE / REAL( MAX( 1, N ) )

         if ( NSIZES.NE.1 ) {
            MTYPES = MIN( MAXTYP, NTYPES )
         } else {
            MTYPES = MIN( MAXTYP+1, NTYPES )
         }

         KA9 = 0
         KB9 = 0
         for (JTYPE = 1; JTYPE <= MTYPES; JTYPE++) { // 640
            IF( .NOT.DOTYPE( JTYPE ) ) GO TO 640
            NMATS = NMATS + 1
            NTEST = 0

            for (J = 1; J <= 4; J++) { // 30
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
               claset('Full', LDA, N, CZERO, CZERO, A, LDA );

            } else if ( ITYPE.EQ.2 ) {

               // Identity

               KA = 0
               KB = 0
               claset('Full', LDA, N, CZERO, CZERO, A, LDA );
               for (JCOL = 1; JCOL <= N; JCOL++) { // 80
                  A( JCOL, JCOL ) = ANORM
   80          CONTINUE

            } else if ( ITYPE.EQ.4 ) {

               // Diagonal Matrix, [Eigen]values Specified

               KA = 0
               KB = 0
               clatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, 0, 0, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE.EQ.5 ) {

               // Hermitian, eigenvalues specified

               KA = MAX( 0, N-1 )
               KB = KA
               clatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, N, N, 'N', A, LDA, WORK, IINFO );

            } else if ( ITYPE.EQ.7 ) {

               // Diagonal, random eigenvalues

               KA = 0
               KB = 0
               clatmr(N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, 0, 0, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE.EQ.8 ) {

               // Hermitian, random eigenvalues

               KA = MAX( 0, N-1 )
               KB = KA
               clatmr(N, N, 'S', ISEED, 'H', WORK, 6, ONE, CONE, 'T', 'N', WORK( N+1 ), 1, ONE, WORK( 2*N+1 ), 1, ONE, 'N', IDUMMA, N, N, ZERO, ANORM, 'NO', A, LDA, IWORK, IINFO );

            } else if ( ITYPE.EQ.9 ) {

               // Hermitian banded, eigenvalues specified

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
               clatms(N, N, 'S', ISEED, 'H', RWORK, IMODE, COND, ANORM, KA, KA, 'N', A, LDA, WORK, IINFO );

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

            // 3) Call CHEGV, CHPGV, CHBGV, CHEGVD, CHPGVD, CHBGVD,
               // CHEGVX, CHPGVX and CHBGVX, do tests.

            // loop over the three generalized problems
                  // IBTYPE = 1: A*x = (lambda)*B*x
                  // IBTYPE = 2: A*B*x = (lambda)*x
                  // IBTYPE = 3: B*A*x = (lambda)*x

            for (IBTYPE = 1; IBTYPE <= 3; IBTYPE++) { // 630

               // loop over the setting UPLO

               for (IBUPLO = 1; IBUPLO <= 2; IBUPLO++) { // 620
                  IF( IBUPLO.EQ.1 ) UPLO = 'U'                   IF( IBUPLO.EQ.2 ) UPLO = 'L'

                  // Generate random well-conditioned positive definite
                  // matrix B, of bandwidth not greater than that of A.

                  clatms(N, N, 'U', ISEED, 'P', RWORK, 5, TEN, ONE, KB, KB, UPLO, B, LDB, WORK( N+1 ), IINFO );

                  // Test CHEGV

                  NTEST = NTEST + 1

                  clacpy(' ', N, N, A, LDA, Z, LDZ );
                  clacpy(UPLO, N, N, B, LDB, BB, LDB );

                  chegv(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, RWORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CHEGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  csgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  // Test CHEGV_2STAGE

                  NTEST = NTEST + 1

                  clacpy(' ', N, N, A, LDA, Z, LDZ );
                  clacpy(UPLO, N, N, B, LDB, BB, LDB );

                  chegv_2stage(IBTYPE, 'N', UPLO, N, Z, LDZ, BB, LDB, D2, WORK, NWORK, RWORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 ) 'CHEGV_2STAGE(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                   // CALL CSGT01( IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z,
      // $                         LDZ, D, WORK, RWORK, RESULT( NTEST ) )

                  // Do Tests | D1 - D2 | / ( |D1| ulp )
                  // D1 computed using the standard 1-stage reduction as reference
                  // D2 computed using the 2-stage reduction

                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  for (J = 1; J <= N; J++) { // 151
                     TEMP1 = MAX( TEMP1, ABS( D( J ) ),  ABS( D2( J ) ) )
                     TEMP2 = MAX( TEMP2, ABS( D( J )-D2( J ) ) )
  151             CONTINUE

                  RESULT( NTEST ) = TEMP2 /  MAX( UNFL, ULP*MAX( TEMP1, TEMP2 ) )

                  // Test CHEGVD

                  NTEST = NTEST + 1

                  clacpy(' ', N, N, A, LDA, Z, LDZ );
                  clacpy(UPLO, N, N, B, LDB, BB, LDB );

                  chegvd(IBTYPE, 'V', UPLO, N, Z, LDZ, BB, LDB, D, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CHEGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  csgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  // Test CHEGVX

                  NTEST = NTEST + 1

                  clacpy(' ', N, N, A, LDA, AB, LDA );
                  clacpy(UPLO, N, N, B, LDB, BB, LDB );

                  chegvx(IBTYPE, 'V', 'A', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CHEGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  csgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1

                  clacpy(' ', N, N, A, LDA, AB, LDA );
                  clacpy(UPLO, N, N, B, LDB, BB, LDB );

                  // since we do not know the exact eigenvalues of this
                  // eigenpair, we just set VL and VU as constants.
                  // It is quite possible that there are no eigenvalues
                  // in this interval.

                  VL = ZERO
                  VU = ANORM
                  chegvx(IBTYPE, 'V', 'V', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CHEGVX(V,V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  csgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1

                  clacpy(' ', N, N, A, LDA, AB, LDA );
                  clacpy(UPLO, N, N, B, LDB, BB, LDB );

                  chegvx(IBTYPE, 'V', 'I', UPLO, N, AB, LDA, BB, LDB, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, NWORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CHEGVX(V,I,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 100
                     }
                  }

                  // Do Test

                  csgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

  100             CONTINUE

                  // Test CHPGV

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 120
                        for (I = 1; I <= J; I++) { // 110
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  110                   CONTINUE
  120                CONTINUE
                  } else {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 140
                        for (I = J; I <= N; I++) { // 130
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  130                   CONTINUE
  140                CONTINUE
                  }

                  chpgv(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, RWORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CHPGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  csgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  // Test CHPGVD

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 160
                        for (I = 1; I <= J; I++) { // 150
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  150                   CONTINUE
  160                CONTINUE
                  } else {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 180
                        for (I = J; I <= N; I++) { // 170
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  170                   CONTINUE
  180                CONTINUE
                  }

                  chpgvd(IBTYPE, 'V', UPLO, N, AP, BP, D, Z, LDZ, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CHPGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  csgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  // Test CHPGVX

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 200
                        for (I = 1; I <= J; I++) { // 190
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  190                   CONTINUE
  200                CONTINUE
                  } else {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 220
                        for (I = J; I <= N; I++) { // 210
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  210                   CONTINUE
  220                CONTINUE
                  }

                  chpgvx(IBTYPE, 'V', 'A', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CHPGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  csgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 240
                        for (I = 1; I <= J; I++) { // 230
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  230                   CONTINUE
  240                CONTINUE
                  } else {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 260
                        for (I = J; I <= N; I++) { // 250
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  250                   CONTINUE
  260                CONTINUE
                  }

                  VL = ZERO
                  VU = ANORM
                  chpgvx(IBTYPE, 'V', 'V', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CHPGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  csgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  NTEST = NTEST + 1

                  // Copy the matrices into packed storage.

                  if ( LSAME( UPLO, 'U' ) ) {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 280
                        for (I = 1; I <= J; I++) { // 270
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  270                   CONTINUE
  280                CONTINUE
                  } else {
                     IJ = 1
                     for (J = 1; J <= N; J++) { // 300
                        for (I = J; I <= N; I++) { // 290
                           AP( IJ ) = A( I, J )
                           BP( IJ ) = B( I, J )
                           IJ = IJ + 1
  290                   CONTINUE
  300                CONTINUE
                  }

                  chpgvx(IBTYPE, 'V', 'I', UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, INFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUNIT, FMT = 9999 )'CHPGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                     INFO = ABS( IINFO )
                     if ( IINFO.LT.0 ) {
                        RETURN
                     } else {
                        RESULT( NTEST ) = ULPINV
                        GO TO 310
                     }
                  }

                  // Do Test

                  csgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

  310             CONTINUE

                  if ( IBTYPE.EQ.1 ) {

                     // TEST CHBGV

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 340
                           DO 320 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  320                      CONTINUE
                           DO 330 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  330                      CONTINUE
  340                   CONTINUE
                     } else {
                        for (J = 1; J <= N; J++) { // 370
                           DO 350 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  350                      CONTINUE
                           DO 360 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  360                      CONTINUE
  370                   CONTINUE
                     }

                     chbgv('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, RWORK, IINFO );
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'CHBGV(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     csgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                     // TEST CHBGVD

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 400
                           DO 380 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  380                      CONTINUE
                           DO 390 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  390                      CONTINUE
  400                   CONTINUE
                     } else {
                        for (J = 1; J <= N; J++) { // 430
                           DO 410 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  410                      CONTINUE
                           DO 420 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  420                      CONTINUE
  430                   CONTINUE
                     }

                     chbgvd('V', UPLO, N, KA, KB, AB, LDA, BB, LDB, D, Z, LDZ, WORK, NWORK, RWORK, LRWORK, IWORK, LIWORK, IINFO );
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'CHBGVD(V,' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     csgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                     // Test CHBGVX

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 460
                           DO 440 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  440                      CONTINUE
                           DO 450 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  450                      CONTINUE
  460                   CONTINUE
                     } else {
                        for (J = 1; J <= N; J++) { // 490
                           DO 470 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  470                      CONTINUE
                           DO 480 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  480                      CONTINUE
  490                   CONTINUE
                     }

                     chbgvx('V', 'A', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'CHBGVX(V,A' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     csgt01(IBTYPE, UPLO, N, N, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 520
                           DO 500 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  500                      CONTINUE
                           DO 510 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  510                      CONTINUE
  520                   CONTINUE
                     } else {
                        for (J = 1; J <= N; J++) { // 550
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
                     chbgvx('V', 'V', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'CHBGVX(V,V' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     csgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                     NTEST = NTEST + 1

                     // Copy the matrices into band storage.

                     if ( LSAME( UPLO, 'U' ) ) {
                        for (J = 1; J <= N; J++) { // 580
                           DO 560 I = MAX( 1, J-KA ), J
                              AB( KA+1+I-J, J ) = A( I, J )
  560                      CONTINUE
                           DO 570 I = MAX( 1, J-KB ), J
                              BB( KB+1+I-J, J ) = B( I, J )
  570                      CONTINUE
  580                   CONTINUE
                     } else {
                        for (J = 1; J <= N; J++) { // 610
                           DO 590 I = J, MIN( N, J+KA )
                              AB( 1+I-J, J ) = A( I, J )
  590                      CONTINUE
                           DO 600 I = J, MIN( N, J+KB )
                              BB( 1+I-J, J ) = B( I, J )
  600                      CONTINUE
  610                   CONTINUE
                     }

                     chbgvx('V', 'I', UPLO, N, KA, KB, AB, LDA, BB, LDB, BP, MAX( 1, N ), VL, VU, IL, IU, ABSTOL, M, D, Z, LDZ, WORK, RWORK, IWORK( N+1 ), IWORK, IINFO );
                     if ( IINFO.NE.0 ) {
                        WRITE( NOUNIT, FMT = 9999 )'CHBGVX(V,I' // UPLO // ')', IINFO, N, JTYPE, IOLDSD
                        INFO = ABS( IINFO )
                        if ( IINFO.LT.0 ) {
                           RETURN
                        } else {
                           RESULT( NTEST ) = ULPINV
                           GO TO 620
                        }
                     }

                     // Do Test

                     csgt01(IBTYPE, UPLO, N, M, A, LDA, B, LDB, Z, LDZ, D, WORK, RWORK, RESULT( NTEST ) );

                  }

  620          CONTINUE
  630       CONTINUE

            // End of Loop -- Check for RESULT(j) > THRESH

            NTESTT = NTESTT + NTEST
            slafts('CSG', N, N, JTYPE, NTEST, RESULT, IOLDSD, THRESH, NOUNIT, NERRS );
  640    CONTINUE
  650 CONTINUE

      // Summary

      slasum('CSG', NOUNIT, NERRS, NTESTT );

      RETURN

 9999 FORMAT( ' CDRVSG2STG: ', A, ' returned INFO=', I6, '.', / 9X, 'N=', I6, ', JTYPE=', I6, ', ISEED=(', 3( I5, ',' ), I5, ')' )

      // End of CDRVSG2STG

      }
