## MEMO
### 命名規則
このプロジェクトは以下のような規則で命名することにする:

|対象|規則|例|
|:--|:--|:--|
|namespace|`lower_camel_case`|`rainbow`|
|class|`UpperCamelCase`|`CudaContext`|
|struct|`UpperCamelCase`|``|
|enum / enum class|`UpperCamelCase`|``|
|型 alias|`UpperCamelCase`|``|
|関数|`lower_snake_case`|``|
|メンバ関数|`lower_snake_case`|``|
|CUDA kernel|`lower_snake_case`|``|
|変数|`lower_snake_case`|``|
|関数の引数|`lower_snake_case`|``|
|public struct field|`lower_snake_case`|``|
|メンバ変数|`lower_snake_case_`|``|
|マクロ|`RAINBOW_UPPER_SNAKE_CASE`|``|
|ファイル名|`lower_snake_case`|``|
|ディレクトリ名|`lower_snake_case`|``|
|CMake target|`rainbow_lower_snake_case`|``|

## クラス，構造体
所有権や不変条件をもつ型は class, そうでなければ struct として実装する