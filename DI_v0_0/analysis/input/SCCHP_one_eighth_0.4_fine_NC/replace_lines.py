#!/usr/bin/env python3
"""input_loc.txt の 7745〜40404 行目を、
second_inverse_new_control_points.txt の 2〜32662 行目で置き換える。
（行番号はすべて 1 始まり・両端を含む）
"""

TARGET_FILE = "input_loc.txt"
SOURCE_FILE = "second_inverse_new_control_points.txt"
OUTPUT_FILE = "input_loc.txt"  # 上書きしたい場合は TARGET_FILE と同じにする

TARGET_START, TARGET_END = 7745, 40404   # 置き換えられる範囲
SOURCE_START, SOURCE_END = 3, 32662      # 置き換えに使う範囲


def main():
    with open(TARGET_FILE, "r", encoding="utf-8", newline="") as f:
        target = f.readlines()
    with open(SOURCE_FILE, "r", encoding="utf-8", newline="") as f:
        source = f.readlines()

    if TARGET_END > len(target):
        raise ValueError(f"{TARGET_FILE} は {len(target)} 行しかありません")
    if SOURCE_END > len(source):
        raise ValueError(f"{SOURCE_FILE} は {len(source)} 行しかありません")

    new_block = source[SOURCE_START - 1:SOURCE_END]

    # 置換範囲の最終行に改行がない場合に備えて、末尾の改行を補う
    if new_block and not new_block[-1].endswith("\n") and TARGET_END < len(target):
        new_block[-1] += "\n"

    result = target[:TARGET_START - 1] + new_block + target[TARGET_END:]

    with open(OUTPUT_FILE, "w", encoding="utf-8", newline="") as f:
        f.writelines(result)

    print(f"置き換え前: {len(target)} 行 "
          f"(うち {TARGET_END - TARGET_START + 1} 行を削除)")
    print(f"挿入: {len(new_block)} 行")
    print(f"置き換え後: {len(result)} 行 -> {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
