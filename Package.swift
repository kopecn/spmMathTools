// swift-tools-version: 6.1
// The swift-tools-version declares the minimum version of Swift required to build this package.

import PackageDescription

let package = Package(
    name: "spmMathTools",
    products: [
        .library(
            name: "spmMathTools",
            targets: ["spmMathTools"]),
    ],
    dependencies: [
        .package(url: "https://github.com/keyvariable/kvSIMD.swift.git", from: "1.1.0"),
        .package(url: "https://github.com/daikimat/depermaid.git", from: "1.1.0"),
    ],
    targets: [
        .target(
            name: "spmMathTools",
            dependencies: [
                .product(name: "kvSIMD", package: "kvSIMD.swift"),
            ],
            path: "spm/Sources/spmMathTools"
        ),
        .testTarget(
            name: "spmMathToolsTests",
            dependencies: ["spmMathTools"],
            path: "spm/Tests/spmMathToolsTests"
        ),
    ]
)
