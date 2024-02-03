      SUBROUTINE ZLAHILB( N, NRHS, A, LDA, X, LDX, B, LDB, WORK, INFO, PATH)

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int     N, NRHS, LDA, LDX, LDB, INFO;
      // .. Array Arguments ..
      double           WORK(N);
      COMPLEX*16 A(LDA,N), X(LDX, NRHS), B(LDB, NRHS)
      String      PATH;
      // ..

*  =====================================================================
      // .. Local Scalars ..
      int     TM, TI, R;
      int     M;
      int     I, J;
      COMPLEX*16 TMP
      String      C2;
      // ..
      // .. Parameters ..
      // NMAX_EXACT   the largest dimension where the generated data is
                   // exact.
      // NMAX_APPROX  the largest dimension where the generated data has
                   // a small componentwise relative error.
      // ??? complex uses how many bits ???
      int     NMAX_EXACT, NMAX_APPROX, SIZE_D;
      const     NMAX_EXACT = 6, NMAX_APPROX = 11, SIZE_D = 8;

      // D's are generated from random permutation of those eight elements.
      COMPLEX*16 D1(8), D2(8), INVD1(8), INVD2(8)
      DATA D1 /(-1.0D0,0.0D0),(0.0D0,1.0D0),(-1.0D0,-1.0D0), (0.0D0,-1.0D0),(1.0D0,0.0D0),(-1.0D0,1.0D0),(1.0D0,1.0D0), (1.0D0,-1.0D0)/       DATA D2 /(-1.0D0,0.0D0),(0.0D0,-1.0D0),(-1.0D0,1.0D0), (0.0D0,1.0D0),(1.0D0,0.0D0),(-1.0D0,-1.0D0),(1.0D0,-1.0D0), (1.0D0,1.0D0)/        DATA INVD1 /(-1.0D0,0.0D0),(0.0D0,-1.0D0),(-0.5D0,0.5D0), (0.0D0,1.0D0),(1.0D0,0.0D0),(-0.5D0,-0.5D0),(0.5D0,-0.5D0), (0.5D0,0.5D0)/       DATA INVD2 /(-1.0D0,0.0D0),(0.0D0,1.0D0),(-0.5D0,-0.5D0), (0.0D0,-1.0D0),(1.0D0,0.0D0),(-0.5D0,0.5D0),(0.5D0,0.5D0), (0.5D0,-0.5D0)/
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. External Functions
      // EXTERNAL ZLASET, LSAMEN
      // INTRINSIC DBLE
      bool    LSAMEN;
      // ..
      // .. Executable Statements ..
      C2 = PATH( 2: 3 )

      // Test the input arguments

      INFO = 0
      if (N .LT. 0 .OR. N .GT. NMAX_APPROX) {
         INFO = -1
      } else if (NRHS .LT. 0) {
         INFO = -2
      } else if (LDA .LT. N) {
         INFO = -4
      } else if (LDX .LT. N) {
         INFO = -6
      } else if (LDB .LT. N) {
         INFO = -8
      }
      if (INFO .LT. 0) {
         CALL XERBLA('ZLAHILB', -INFO)
         RETURN
      }
      if (N .GT. NMAX_EXACT) {
         INFO = 1
      }

      // Compute M = the LCM of the integers [1, 2*N-1].  The largest
      // reasonable N is small enough that integers suffice (up to N = 11).
      M = 1
      DO I = 2, (2*N-1)
         TM = M
         TI = I
         R = MOD(TM, TI)
         DO WHILE (R .NE. 0)
            TM = TI
            TI = R
            R = MOD(TM, TI)
         END DO
         M = (M / TI) * I
      END DO

      // Generate the scaled Hilbert matrix in A
      // If we are testing SY routines,
        t // ake D1_i = D2_i, else, D1_i = D2_i*
      if ( LSAMEN( 2, C2, 'SY' ) ) {
         DO J = 1, N
            DO I = 1, N
               A(I, J) = D1(MOD(J,SIZE_D)+1) * (DBLE(M) / (I + J - 1)) * D1(MOD(I,SIZE_D)+1)
            END DO
         END DO
      } else {
         DO J = 1, N
            DO I = 1, N
               A(I, J) = D1(MOD(J,SIZE_D)+1) * (DBLE(M) / (I + J - 1)) * D2(MOD(I,SIZE_D)+1)
            END DO
         END DO
      }

      // Generate matrix B as simply the first NRHS columns of M * the
      // identity.
      TMP = DBLE(M)
      CALL ZLASET('Full', N, NRHS, (0.0D+0,0.0D+0), TMP, B, LDB)

      // Generate the true solutions in X.  Because B = the first NRHS
      // columns of M*I, the true solutions are just the first NRHS columns
      // of the inverse Hilbert matrix.
      WORK(1) = N
      DO J = 2, N
         WORK(J) = (  ( (WORK(J-1)/(J-1)) * (J-1 - N) ) /(J-1)  ) * (N +J -1)
      END DO

      // If we are testing SY routines,
           t // ake D1_i = D2_i, else, D1_i = D2_i*
      if ( LSAMEN( 2, C2, 'SY' ) ) {
         DO J = 1, NRHS
            DO I = 1, N
               X(I, J) = INVD1(MOD(J,SIZE_D)+1) * ((WORK(I)*WORK(J)) / (I + J - 1)) * INVD1(MOD(I,SIZE_D)+1)
            END DO
         END DO
      } else {
         DO J = 1, NRHS
            DO I = 1, N
               X(I, J) = INVD2(MOD(J,SIZE_D)+1) * ((WORK(I)*WORK(J)) / (I + J - 1)) * INVD1(MOD(I,SIZE_D)+1)
            END DO
         END DO
      }
      }
