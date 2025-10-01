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
    targets: [
        .target(
            name: "spmMathTools",
            path: "spm/Sources/spmMathTools"
        ),
        .testTarget(
            name: "spmMathToolsTests",
            dependencies: ["spmMathTools"],
            path: "spm/Tests/spmMathToolsTests"
        ),
    ]
)
