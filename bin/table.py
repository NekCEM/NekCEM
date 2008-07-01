class Table:
    def __init__(self):
        self.Rows = []

    def append_row(self, entries):
        self.Rows.append(entries)

    def stringify(self, val, row, col):
        return str(val)

    def entry(self, row, col):
        if col >= len(self.Rows[row]):
            return ""
        else:
            return self.stringify(self.Rows[row][col], row, col)

    def __str__(self):
        num_columns = max([len(i) for i in self.Rows])
        column_widths = [
          max([len(self.entry(row, col))
                  for row in range(len(self.Rows))])
          for col in range(num_columns)]

        def filled_entry(row, col):
            result = self.entry(row, col)
            result += (column_widths[col]-len(result))*" "
            return result

        def format_row(row):
            return "|".join([filled_entry(row, col) 
                    for col in range(num_columns)])

        separator = "+".join([column_widths[col]*"-"
                    for col in range(num_columns)])

        return "\n".join([format_row(0), separator] +
                [format_row(i) for i in range(1, len(self.Rows))])



if __name__ == "__main__":
    t = Table()
    t.append_row(["yes", "no"])
    t.appenddd_row(["yes", "no", "hullo"])
    t.appenddd_row(["yessa", "no", "hullo"])
    t.appenddd_row(["yessa", "nona", "hullo"])
    print t
