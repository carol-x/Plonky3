// credit to https://github.com/HorizenLabs/poseidon2/blob/main/plain_implementations/src/poseidon2/poseidon2.rs
use super::poseidon2_params::Poseidon2Params;
use crate::fields::goldilocks::FpGoldiLocks;
use crate::fields::utils::from_hex;

use lazy_static::lazy_static;
use std::sync::Arc;

type Scalar = FpGoldiLocks;

lazy_static! {
    pub static ref MAT_DIAG8_M_1: Vec<Scalar> = vec![
    from_hex("0xd57b33d215cc4805"),
    from_hex("0xaa2238eb3ac17b62"),
    from_hex("0x28925fe2f3895c0d"),
    from_hex("0x3dab9370a67db22e"),
    from_hex("0xe5cafe41ef4eac62"),
    from_hex("0x4c633d43f2260c06"),
    from_hex("0x1fa5fb8a31d6369d"),
    from_hex("0x999a460e4a706453"),
    ];

    pub static ref MAT_INTERNAL8: Vec<Vec<Scalar>> = vec![
    vec![from_hex("0xd57b33d215cc4806"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    ],
    vec![from_hex("0x0000000000000001"),
    from_hex("0xaa2238eb3ac17b63"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    ],
    vec![from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x28925fe2f3895c0e"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    ],
    vec![from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x3dab9370a67db22f"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    ],
    vec![from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0xe5cafe41ef4eac63"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    ],
    vec![from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x4c633d43f2260c07"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    ],
    vec![from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x1fa5fb8a31d6369e"),
    from_hex("0x0000000000000001"),
    ],
    vec![from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x0000000000000001"),
    from_hex("0x999a460e4a706454"),
    ],
    ];

    pub static ref RC8: Vec<Vec<Scalar>> = vec![
    vec![from_hex("0x57056152cedf0fe7"),
    from_hex("0x44b125d16e93ca85"),
    from_hex("0x8e8ea2ff8b7a6d2a"),
    from_hex("0xcce7c6cc1468fa13"),
    from_hex("0x47f5feb953ce5073"),
    from_hex("0xfd8f41d8ee6b700e"),
    from_hex("0xe40f59b8db57aeb7"),
    from_hex("0x78b572234ff68244"),
    ],
    vec![from_hex("0x926b547a9712ed0b"),
    from_hex("0xb1525da069ba226c"),
    from_hex("0xf37650e9d8ef46d3"),
    from_hex("0x3146518c7738aefc"),
    from_hex("0x04aa9f4d916e9e5b"),
    from_hex("0xde603b81bb63d21c"),
    from_hex("0x8382c29e88cf2c81"),
    from_hex("0x50456f59f404cb88"),
    ],
    vec![from_hex("0x44bda4a6711f6ddb"),
    from_hex("0xe4c94cbc9e7d15b7"),
    from_hex("0x7faec52ce37a8256"),
    from_hex("0x7748e71fd7803107"),
    from_hex("0x9b6baf83e49be593"),
    from_hex("0xd47fe8a5c8b27ed3"),
    from_hex("0xfcdf1e28d16392ad"),
    from_hex("0x976753b4b516a9ee"),
    ],
    vec![from_hex("0xc16ea705aa7ee467"),
    from_hex("0x18183d87f912ebbb"),
    from_hex("0x02d3b175b21777fe"),
    from_hex("0x98e4c2d93e0aaaef"),
    from_hex("0xc31191d90cd41c96"),
    from_hex("0x69f8f94595ad453e"),
    from_hex("0x1de4127f3e248a2d"),
    from_hex("0xbcce9849c99a069c"),
    ],
    vec![from_hex("0x8b8e707932590779"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x4d7fff707c77890f"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x7d36116962851777"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x1dc9f40fbb3146b7"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x6a235e2d5bef54e0"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x4d1a9ae6dd337207"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x46ab49a6009cda1a"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x78e759e819648587"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0xee6e84b7763598a4"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x0b426bdcaad3050e"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x1f3cd981be91490e"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0xd54572f7ecf947a1"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x393c4432d0e86a1e"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x3f1b43149ef3f4f8"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x3705f6a66d25dce4"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x3e809302b3d41471"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x6e50830e082b17f1"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x711232bf2d77ac38"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x4235f7d079c78096"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0xab1bbdc696a72a25"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0xdb1ef6f3f7fed243"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0xd21981014e77d809"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    from_hex("0x0000000000000000"),
    ],
    vec![from_hex("0x5b2cb2bd03a18856"),
    from_hex("0x8e45a3e4bf30df6c"),
    from_hex("0x3f9948080379716d"),
    from_hex("0x41c2ba50c09d6c70"),
    from_hex("0x5c2f57c6f81d2c6b"),
    from_hex("0x91cfb3d3b4b04a7a"),
    from_hex("0x81327090650355f6"),
    from_hex("0x06957eabf4817942"),
    ],
    vec![from_hex("0x7f08201e9da0e064"),
    from_hex("0x7467dfc268e1d6e0"),
    from_hex("0x38a9992ed589cc80"),
    from_hex("0x266a6e035fee9286"),
    from_hex("0xd19ebfbf75ffbf79"),
    from_hex("0x9f1dc0303ca0acfb"),
    from_hex("0x230f2d6a36b23347"),
    from_hex("0xde0cdaab08319a52"),
    ],
    vec![from_hex("0xff9e2984d5f675ba"),
    from_hex("0x27a10c5aca2fcf50"),
    from_hex("0x8982ec2da08deb87"),
    from_hex("0x89f9b8d33e98a684"),
    from_hex("0x269bcee2edb77b24"),
    from_hex("0xcd7fb3f592ab464f"),
    from_hex("0x05060bc8d4341e72"),
    from_hex("0xa75ab333263a6658"),
    ],
    vec![from_hex("0x3962fe1b4bb486e7"),
    from_hex("0x52160689b78a2fd1"),
    from_hex("0x9e953026b7be93e6"),
    from_hex("0x7215465ca2fa2b5a"),
    from_hex("0x458b8385c2107d5b"),
    from_hex("0xd86fd0264024aad9"),
    from_hex("0x2cb61942ee72b44c"),
    from_hex("0x50784c715273f7e7"),
    ],
    ];

    pub static ref POSEIDON2_GOLDILOCKS_8_PARAMS: Arc<Poseidon2Params<Scalar>> = Arc::new(Poseidon2Params::new(8, 7, 8, 22, &MAT_DIAG8_M_1, &MAT_INTERNAL8, &RC8)
    );

}