import numpy as np

def wpgma(
  matrix,
  labels):
    np.fill_diagonal(matrix, np.inf)
    labels = [[label, [label], 0] for label in labels]
    while len(labels) > 1:
        i, j = divmod(
          matrix.argmin(),
          matrix.shape[1])
        minval = matrix[i, j]

        if minval == 0:
            labels[i][0] += " | " + labels[j][0]
        else:
            dist = minval/2
            labels[i][0] = (
              f"({labels[i][0]}:{dist-labels[i][2]},"
              f"{labels[j][0]}:{dist-labels[j][2]})")
            labels[i][2] = dist

            matrix[i] = [
              (matrix[i, col] + matrix[j, col])/2 
              if col != i else np.inf
              for col in range(matrix.shape[1])]
            matrix[:, i] = [
              (matrix[row, i] + matrix[row, j])/2
              if row != i else np.inf
              for row in range(matrix.shape[0])]

        labels[i][1] += labels.pop(j)[1]
        matrix = np.delete(matrix, j, axis=0)
        matrix = np.delete(matrix, j, axis=1)

    return labels[0][0] + ":0.0;", labels[0][1]
