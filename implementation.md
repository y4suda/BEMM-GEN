# CEMM-gen の実装

## TODO
* !githubでの管理の練習
* !ディレクトリ構造の決定
* !処理の流れを考える
* 必要なライブラリを探す
  * 引数の処理
  * setuptools
  * ドキュメント化
* チュートリアル作成

## Github での管理の練習
基本的にはVS codeを使えば良い
1. 左のエクスプローラーからリポジトリの複製をして、新しいウィンドウを作る
1. 適当な作業ディレクトリ内にクローンされる
1. ファイルを編集する
1. 左の「ソース管理」に変更点がリストアップされる
1. 「+」を押してステージに上げる（やめるときは「-」でステージから下ろす）
1. commit文をかいてcommit、同期する
1. 次の作業前に必ず同期すること！

最初にディレクトリ内で設定が必要？

```
git config --global user.name "username"
git config --global user.email "mail@address"
```

## 基本的な処理の流れ

## ディレクトリ構造

## ライブラリ選定

### 引数の処理
標準ライブラリの argparse を使用

サブコマンドを分けてパーサーを用意
https://qiita.com/kzkadc/items/e4fc7bc9c003de1eb6d0
https://qiita.com/oohira/items/308bbd33a77200a35a3d


### インストール方法
https://packaging.python.org/ja/latest/guides/distributing-packages-using-setuptools/
https://qiita.com/propella/items/5cd89caee6379920d889

### ドキュメント化
Sphinxがいい？
https://zenn.dev/yamagishihrd/articles/8b5279a07000a4

docstring のスタイルはGoogle style と Numpy styleがある→とりあえず google style?