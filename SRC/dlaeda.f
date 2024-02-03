      SUBROUTINE DLAEDA( N, TLVLS, CURLVL, CURPBM, PRMPTR, PERM, GIVPTR, GIVCOL, GIVNUM, Q, QPTR, Z, ZTEMP, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                CURLVL, CURPBM, INFO, N, TLVLS;
      // ..
      // .. Array Arguments ..
      int                GIVCOL( 2, * ), GIVPTR( * ), PERM( * ), PRMPTR( * ), QPTR( * );
      double             GIVNUM( 2, * ), Q( * ), Z( * ), ZTEMP( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, HALF, ONE;
      const              ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      int                BSIZ1, BSIZ2, CURR, I, K, MID, PSIZ1, PSIZ2, PTR, ZPTR1;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMV, DROT, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, INT, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0

      if ( N.LT.0 ) {
         INFO = -1
      }
      if ( INFO.NE.0 ) {
         xerbla('DLAEDA', -INFO );
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Determine location of first number in second half.

      MID = N / 2 + 1

      // Gather last/first rows of appropriate eigenblocks into center of Z

      PTR = 1

      // Determine location of lowest level subproblem in the full storage
      // scheme

      CURR = PTR + CURPBM*2**CURLVL + 2**( CURLVL-1 ) - 1

      // Determine size of these matrices.  We add HALF to the value of
      // the SQRT in case the machine underestimates one of these square
      // roots.

      BSIZ1 = INT( HALF+SQRT( DBLE( QPTR( CURR+1 )-QPTR( CURR ) ) ) )
      BSIZ2 = INT( HALF+SQRT( DBLE( QPTR( CURR+2 )-QPTR( CURR+1 ) ) ) )
      DO 10 K = 1, MID - BSIZ1 - 1
         Z( K ) = ZERO
      } // 10
      dcopy(BSIZ1, Q( QPTR( CURR )+BSIZ1-1 ), BSIZ1, Z( MID-BSIZ1 ), 1 );
      dcopy(BSIZ2, Q( QPTR( CURR+1 ) ), BSIZ2, Z( MID ), 1 );
      DO 20 K = MID + BSIZ2, N
         Z( K ) = ZERO
      } // 20

      // Loop through remaining levels 1 -> CURLVL applying the Givens
      // rotations and permutation and then multiplying the center matrices
      // against the current Z.

      PTR = 2**TLVLS + 1
      DO 70 K = 1, CURLVL - 1
         CURR = PTR + CURPBM*2**( CURLVL-K ) + 2**( CURLVL-K-1 ) - 1
         PSIZ1 = PRMPTR( CURR+1 ) - PRMPTR( CURR )
         PSIZ2 = PRMPTR( CURR+2 ) - PRMPTR( CURR+1 )
         ZPTR1 = MID - PSIZ1

        // Apply Givens at CURR and CURR+1

         DO 30 I = GIVPTR( CURR ), GIVPTR( CURR+1 ) - 1
            drot(1, Z( ZPTR1+GIVCOL( 1, I )-1 ), 1, Z( ZPTR1+GIVCOL( 2, I )-1 ), 1, GIVNUM( 1, I ), GIVNUM( 2, I ) );
         } // 30
         DO 40 I = GIVPTR( CURR+1 ), GIVPTR( CURR+2 ) - 1
            drot(1, Z( MID-1+GIVCOL( 1, I ) ), 1, Z( MID-1+GIVCOL( 2, I ) ), 1, GIVNUM( 1, I ), GIVNUM( 2, I ) );
         } // 40
         PSIZ1 = PRMPTR( CURR+1 ) - PRMPTR( CURR )
         PSIZ2 = PRMPTR( CURR+2 ) - PRMPTR( CURR+1 )
         DO 50 I = 0, PSIZ1 - 1
            ZTEMP( I+1 ) = Z( ZPTR1+PERM( PRMPTR( CURR )+I )-1 )
         } // 50
         DO 60 I = 0, PSIZ2 - 1
            ZTEMP( PSIZ1+I+1 ) = Z( MID+PERM( PRMPTR( CURR+1 )+I )-1 )
         } // 60

         // Multiply Blocks at CURR and CURR+1

         // Determine size of these matrices.  We add HALF to the value of
         // the SQRT in case the machine underestimates one of these
         // square roots.

         BSIZ1 = INT( HALF+SQRT( DBLE( QPTR( CURR+1 )-QPTR( CURR ) ) ) )
         BSIZ2 = INT( HALF+SQRT( DBLE( QPTR( CURR+2 )-QPTR( CURR+ 1 ) ) ) )
         if ( BSIZ1.GT.0 ) {
            dgemv('T', BSIZ1, BSIZ1, ONE, Q( QPTR( CURR ) ), BSIZ1, ZTEMP( 1 ), 1, ZERO, Z( ZPTR1 ), 1 );
         }
         dcopy(PSIZ1-BSIZ1, ZTEMP( BSIZ1+1 ), 1, Z( ZPTR1+BSIZ1 ), 1 );
         if ( BSIZ2.GT.0 ) {
            dgemv('T', BSIZ2, BSIZ2, ONE, Q( QPTR( CURR+1 ) ), BSIZ2, ZTEMP( PSIZ1+1 ), 1, ZERO, Z( MID ), 1 );
         }
         dcopy(PSIZ2-BSIZ2, ZTEMP( PSIZ1+BSIZ2+1 ), 1, Z( MID+BSIZ2 ), 1 );

         PTR = PTR + 2**( TLVLS-K )
      } // 70

      RETURN

      // End of DLAEDA

      }
